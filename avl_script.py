import numpy as np
import subprocess
from scipy.optimize import differential_evolution, NonlinearConstraint
import re
import json
import logging
import uuid
import tempfile
import time

# Optimizing with spec_cl - max_efficiency + spec_cl - target_cl leads avl to raise the alpha indefinetly giving an error,
# likely due to not having a limit on the Sref, optimizing without limit on Sref needs more attention

# Useful wrapper that allows to time a function by adding @timer in the line before the function definition
def timer(func):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        print(f"{func.__name__} took {end-start:.4f} seconds")
    return wrapper

"""
TO-DO
possibly add a penalization for bending moment
incorporate polar in .avl file
add regulation for # cores and .json path
add min drag op_type
"""
# Took 1 h 10 min to run
# # 37 min with linear twist
# AN ACCURATE POLAR OF THE PROFILE IS KEY FOR THE OPTIMIZATION

# Calculates Xle for each point in y_pan with a sweep angle of dihed in DEGREES
def CalculateXle(chords:list, y_pan:np.ndarray, sweep=0.0):
    x_pan = np.array([chords[0] - c for c in chords])
    if sweep != 0.0:
        x_pan += y_pan* np.sin(np.pi/180 *sweep)  
    return x_pan

# Calculates Yle for each section with cosine distribution
def CalculateYle(b2:float, N:int):
    y_pan = b2*np.cos(np.pi/2 * np.linspace(0, N-1, N)/(N-1))
    return y_pan[::-1]


# Calculates Zle for each point in y_pan with a dihedral angle of dihed in DEGREES
def CalculateZle(dihed:float, y_pan:np.ndarray):
    z_pan = y_pan * np.sin(np.pi/180 *dihed)
    return z_pan


# Calculates the mean aerodynamic chord
def CalculateMAC(Sref:float, y_pan:np.ndarray, chords:list):
    if len(chords) != len(y_pan):
        raise ValueError("Chords and y_pan don't have the same length")
    chord_squared = np.square(np.array(chords))
    Mac = 2/Sref * np.trapz(chord_squared, y_pan)
    return Mac

# Calculates the wing surface
def CalculateSref(y_pan:np.ndarray, chords:list):
    if len(chords) != len(y_pan):
        raise ValueError("Chords and y_pan don't have the same length")
    Sref =  2 * np.trapz(chords, y_pan)
    return Sref

# Calculates twist for each section given the twist angle in DEGREES
def CalculateTwists(twist:float, y_pan:np.ndarray, b2:float):
    twists = twist/b2 * y_pan
    return twists


def ReadConfigJSON(config_path:str):
    with open(config_path, 'r') as f:
        config = json.load(f)

    file_name    = config["file_name"]
    wing_name    = config["wing_name"]
    profile_name = config["profile_file"]
    CDp          = config["CDp"]
    sections     = config["num_sect"]
    Starget = config["Starget"]
    opt_points = config["opt_points"]

    # [span_b, dihed_b, twist_b, chords_b]
    bounds = []
    bounds.append(config["span"]["bounds"])
    bounds.append(config["dihed"]["bounds"])
    bounds.append(config["twist"]["bounds"])

    # [span, dihed, twist, chords]
    x0 = []
    x0.append(config["span"]["initial"])
    x0.append(config["dihed"]["initial"])
    x0.append(config["twist"]["initial"])

    chords = []
    chords_b = []

    for section in config["sections"]:
        chord_data =section["chord"]
        chords.append(chord_data["initial"])
        chords_b.append(tuple(chord_data["bounds"]))

    x0.extend(chords)
    bounds.extend(chords_b)
    
    return file_name, wing_name, profile_name, CDp, sections, Starget, x0, bounds, opt_points
    

# Composes the .avl file for a surface
def WriteAVLText(
        file_name:str,
        wing_name:str, 
        Ma:float, 
        sections:int,
        chords:list, 
        twists:list, 
        profile_file:str,
        AVL_path:str,
        Xle:np.ndarray, Yle:np.ndarray, Zle:np.ndarray,
        Sref:float, Cref:float, Bref:float, 
        x=0.0, y=0.0, z=0.0,
        Xref=0.0, Yref=0.0, Zref=0.0,
        CDp=0.020
        ):
    BaseAVLText =  f"""{file_name}
{Ma:.2f}                   !   Mach
0     0     0.0       !   iYsym  iZsym  Zsym
{Sref:.4f} {Cref:.4f}  {Bref:.4f}       !   Sref   Cref   Bref   reference area, chord, span
{Xref}  {Yref}   {Zref}       !   Xref   Yref   Zref   moment reference location (arb.)
{CDp}                 !   CDp
#
#==============================================================
#
SURFACE
{wing_name}
6  1.0  {4*sections}  -2.0   ! Nchord   Cspace   Nspan  Sspace
#
# reflect image wing about y=0 plane
YDUPLICATE
     0.00000 
#
# twist angle bias for whole surface
ANGLE
     0.00000    
#
SCALE
  1.0   1.0   1.0
#
# x,y,z bias for whole surface
TRANSLATE
    {x:.10f}    {y:.10f}     {z:.10f}
"""
    for i in range(sections):
        SectionAVLText = f"""#
#--------------------------------------------------------------
#    Xle         Yle         Zle         chord       angle   Nspan  Sspace
SECTION
    {Xle[i]:.10f}         {Yle[i]:.10f}         {Zle[i]:.10f}         {chords[i]:.10f}        {twists[i]:.10f}   4     -2

AFIL
{profile_file}
"""

        BaseAVLText += SectionAVLText

    BaseAVLText += f"""#
#==============================================================
#
"""
    return BaseAVLText


# Runs AVL from the avl_geom_text treated as an .avl file with constraint on alpha to have Cl_target, outputs to out_path
def RunAVL(op_mode:str, target:float, avl_geom_text:str)->str:
    # Create a temporary .avl file
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.avl', delete=False) as temp_geom:
        temp_geom.write(avl_geom_text)
        temp_geom_path = temp_geom.name

    # Define the command sequence you want to send to AVL
    if op_mode == 'spec_cl':
        avl_commands = f"""
load {temp_geom_path}
oper
A
C {target}
x

quit

"""  
    elif op_mode == 'spec_al':
        avl_commands = f"""
load {temp_geom_path}
oper
A
A {target}
x

quit

"""

    # Start the AVL process
    process = subprocess.Popen(
        ['avl'],  
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True  # Ensure strings are used instead of bytes
    )

    # Send the command string and capture output
    stdout, stderr = process.communicate(avl_commands)

    if process.returncode != 0:
        raise RuntimeError(f"AVL failed with return code {process.returncode}")
    
    return stdout


# Gets Cl/Cd from loads_file, use ParseAVLstdout as it is faster
def ParseClCdAVL(loads_file:str):

    with open(loads_file, 'r') as file:
        text = file.read()

    CL_match = re.search(r"CLtot\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)
    CD_match = re.search(r"CDtot\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)

    if not CL_match or not CD_match:
        raise ValueError("CL or CD not found in loads file.")
    
    CLtot = float(CL_match.group(1))
    CDtot = float(CD_match.group(1))

    return CLtot, CDtot


# Gets Cl/Cd from AVL's stdout
def ParseAVLstdout(stdout:str):
    cltot = cdtot = None
    cl_pattern = re.compile(r'^\s*CLtot\s*=\s*([-+]?[\d.]+)')
    cd_pattern = re.compile(r'^\s*CDtot\s*=\s*([-+]?[\d.]+)')

    for line in stdout.splitlines():
        if cltot is None:
            cl_match = cl_pattern.match(line)
            if cl_match:
                cltot = float(cl_match.group(1))
        if cdtot is None:
            cd_match = cd_pattern.match(line)
            if cd_match:
                cdtot = float(cd_match.group(1))
        if cltot is not None and cdtot is not None:
            break

    if cltot is None:
        raise RuntimeError('AVL output missing CLtot')
    if cdtot is None:
        raise RuntimeError('AVL output missing CDtot')

    return cltot, cdtot

def RunAVLSimulations(params:list, fixed_params:dict):

    uid = uuid.uuid4().hex  # generate unique ID for each call
    loads_file = f"temp/loads_{uid}.txt"
    AVL_path   = f"temp/input_{uid}.avl"
    
    Bref         = params[0]
    dihed        = params[1]
    twist        = params[2]
    chords       = params[3:]
    
    sections     = fixed_params['sections']
    file_name    = fixed_params['file_name']
    wing_name    = fixed_params['wing_name']
    profile_file = fixed_params['profile_file']
    CDp          = fixed_params['CDp']
    opt_points   = fixed_params['opt_points']


    Yle = CalculateYle(Bref/2, sections)
    Xle = CalculateXle(chords, Yle)
    Zle = CalculateZle(dihed, Yle)

    twists = CalculateTwists(twist, Yle, Bref/2)
    
    Sref = CalculateSref(Yle, chords)
    Cref = CalculateMAC(Sref, Yle, chords)

    LD_ratio = 0.0
    CLs = []

    for point in opt_points:
        mode = point['op_mode']
        op_type = point['op_type']
        
        Ma = point['Ma']
        target = (point['op_point'] / Sref) if (mode == 'spec_cl') else point['op_point']

        avl_text = WriteAVLText(
            file_name, 
            wing_name, 
            Ma, 
            sections, 
            chords, 
            twists, 
            profile_file, 
            AVL_path, 
            Xle, Yle, Zle,
            Sref, Cref, Bref,
            0, 0, 0,
            0, 0, 0,
            CDp)
        
        if mode == 'spec_cl' and op_type == 'target_cl':
                raise ValueError("op_mode 'spec_cl' is not compatible with op_type 'target_cl'")
            
        stdout = RunAVL(mode, target, avl_text)
        cl, cd = ParseAVLstdout(stdout)

        if op_type == 'max_efficiency':
            LD_ratio += point['weight'] * (cl/cd)
        elif op_type == 'target_cl':
            CLs.append(cl)

    logging.info(f'\nCl/Cd: {LD_ratio}')

    return -LD_ratio, CLs


class AVLEvaluator():
    def __init__(self):
        self.cache = {}
    
    def evaluate(self, params, fixed_params):
        results = RunAVLSimulations(params, fixed_params)
        self.cache['last'] = results
        return results[0]
    
    def get_CLS(self):
        return self.cache['last'][1]

      
def WriteFinalResults(params:list, fixed_params:dict):

    Bref         = params[0]
    dihed        = params[1]
    twist        = params[2]
    chords       = params[3:]
    
    sections     = fixed_params['sections']
    file_name    = fixed_params['file_name']
    wing_name    = fixed_params['wing_name']
    profile_file = fixed_params['profile_file']
    CDp          = fixed_params['CDp']

    opt_points   = fixed_params['opt_points']
    Ma = opt_points[0]['Ma']

    AVL_path     = 'Results.avl'

    Yle = CalculateYle(Bref/2, sections)
    Xle = CalculateXle(chords, Yle)
    Zle = CalculateZle(dihed, Yle)

    twists = CalculateTwists(twist, Yle, Bref/2)

    Sref = CalculateSref(Yle, chords)
    Cref = CalculateMAC(Sref, Yle, chords)

    final_text = WriteAVLText(
        file_name,
        wing_name,
        Ma,
        sections,
        chords,
        twists,
        profile_file,
        AVL_path,
        Xle, Yle, Zle,
        Sref, Cref, Bref,
        CDp=CDp
    )
    with open('Results.avl', 'w') as r:
        r.write(final_text)

# Using global variable to make cost function picklable
evaluator = AVLEvaluator()
def CostFunction(params:list, fixed_params:dict):
        return evaluator.evaluate(params, fixed_params)


def main():
    config_path = "config.json"
    
    (
        file_name, wing_name, profile_name, 
        CDp, sections, Starget,
        x0, bounds, opt_points
    ) = ReadConfigJSON(config_path)

    fixed_params = {
        'file_name': file_name,
        'wing_name': wing_name,
        'profile_file': profile_name,
        'CDp': CDp,
        'sections': sections,
        'opt_points':opt_points
    }

    cl_targets = np.array([opt_point['target_cl'] for opt_point in opt_points if opt_point['op_mode'] == 'spec_al'])
    cl_tol = 0.05

    evaluator.evaluate(x0, fixed_params)
    
    def CLsConstraint(params:list):
        return evaluator.get_CLS()
    
    def SrefConstraint(params:list):
        Yle = CalculateYle(params[0]/2, sections)
        S_true = CalculateSref(Yle, params[3:])
        return S_true - Starget
    
    cl_constraints = NonlinearConstraint(CLsConstraint, cl_targets - cl_tol, cl_targets + cl_tol)
    sref_constraint = NonlinearConstraint(SrefConstraint,0.5*Starget, 1.5*Starget )

    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')
    
    res = differential_evolution(
        CostFunction,
        bounds=bounds,
        args=(fixed_params,),
        strategy="best1bin",
        popsize=30,                 # Number of candidates for generation
        atol=0.00001,
        maxiter=100,                # Number of generations, 100 gave good results
        polish=False,               # Uses gradient based method for final local search, gives some improvement but takes forever (hours)
        workers=-1,                 # Adjust to number of cores you want to use, -1 for all aviable
        x0 = x0,
        disp=True,
        updating='deferred',
        mutation=(0.5, 1.5),
        recombination=0.7, 
        constraints=[cl_constraints, sref_constraint]
    )

    logging.info(res)

    WriteFinalResults(res.x, fixed_params)


if __name__ == '__main__':
    main()
