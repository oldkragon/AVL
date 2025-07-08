import numpy as np
import subprocess
from scipy.optimize import differential_evolution, NonlinearConstraint
import re
import json
import logging
import uuid

"""
TO-DO
insert a mode for linear twist and test if gives better results careful with ReadConfigJSON, CostFunction, WriteFinalResults
still needs some work
multiprocessing doesn't read file correctly (same S for each run in iteration)
add regulation for # cores and .json path
"""
# Took 1 h 10 min to run, see with linear twist
# AN ACCURATE ESTIMATE OF THE PROFILE CD AT 0-1 deg IS KEY FOR THE OPTIMIZATION

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

    Starget      = config["Starget"]
    Bref         = config["Bref"]
    Ma           = config["Ma"]
    CLtarget     = config["CLtarget"]
    CDp          = config["CDp"]

    sections = config["num_sect"]

    # [dihed_b, chords_b, twists_b]
    bounds = []
    # [dihed, chords, twists]
    x0 = []

    chords = []
    twists = []

    chords_b = []
    twists_b = []

    dihed_data = config["dihed"]
    bounds.append(dihed_data["bounds"])

    dihed = dihed_data["initial"]
    x0.append(dihed)

    for section in config["sections"]:
        chord_data =section["chord"]
        chords.append(chord_data["initial"])
        chords_b.append(tuple(chord_data["bounds"]))

        twist_data =section["twist"]
        twists.append(twist_data["initial"])
        twists_b.append(tuple(twist_data["bounds"]))

    x0.extend(chords)
    x0.extend(twists)

    bounds.extend(chords_b)
    bounds.extend(twists_b)

    bounds = bounds
    
    return file_name, wing_name, profile_name, Starget, Bref, Ma, CLtarget, CDp, sections, x0, bounds
    

# Composes the .avl file for a surface
def WriteAVLFile(
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
10  1.0  26  1.0   ! Nchord   Cspace   Nspan  Sspace
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
    with open(AVL_path, 'w') as avl_file:
        avl_file.write(BaseAVLText)
    return BaseAVLText


# Runs AVL from the in_path .avl file with constraint on alpha to have Cl_target, outputs to out_path
def RunAVL(Cl_target:float, AVL_path:str, out_path:str):

    # Define the command sequence you want to send to AVL
    avl_commands = f"""
load {AVL_path}
oper
A
C {Cl_target}
x
w
{out_path}

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

    logging.debug("AVL STDOUT:\n" + stdout)
    logging.debug("AVL STDERR:\n" + stderr)

    if process.returncode != 0:
        raise RuntimeError(f"AVL failed with return code {process.returncode}")



# Gets Cl/Cd from loads_file
def ParseClCd(loads_file:str):

    with open(loads_file, 'r') as file:
        text = file.read()

    CL_match = re.search(r"CLtot\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)
    CD_match = re.search(r"CDtot\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)

    if not CL_match or not CD_match:
        raise ValueError("CL or CD not found in loads file.")
    
    CLtot = float(CL_match.group(1))
    CDtot = float(CD_match.group(1))

    return CLtot/CDtot


def CostFunction(params:list, fixed_params:dict):

    uid = uuid.uuid4().hex  # generate unique ID for each call
    loads_file = f"temp/loads_{uid}.txt"
    AVL_path   = f"temp/input_{uid}.avl"
    
    dihed        = params[0]
    sections     = fixed_params['sections']
    chords       = params[1:1 + sections]
    twists       = params[1 + sections:]

    file_name    = fixed_params['file_name']
    wing_name    = fixed_params['wing_name']
    profile_file = fixed_params['profile_file']
    Bref         = fixed_params['Bref']
    Ma           = fixed_params['Ma']
    CLtarget     = fixed_params['CLtarget']
    Starget      = fixed_params['Starget']
    CDp          = fixed_params['CDp']


    Yle = CalculateYle(Bref/2, sections)
    Xle = CalculateXle(chords, Yle)
    Zle = CalculateZle(dihed, Yle)
    
    Sref = CalculateSref(Yle, chords)
    Cref = CalculateMAC(Sref, Yle, chords)

    WriteAVLFile(
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
    
    RunAVL(CLtarget, AVL_path, loads_file)

    LD_ratio = ParseClCd(loads_file)

    logging.info(f"\nParams: {params}\n -L/D: {-LD_ratio}\n Sref : {Sref}")

    return -LD_ratio


def WriteFinalResults(params:list, fixed_params:dict):
    dihed        = params[0]
    sections     = fixed_params['sections']
    chords       = params[1:1 + sections]
    twists       = params[1 + sections:]

    file_name    = fixed_params['file_name']
    wing_name    = fixed_params['wing_name']
    profile_file = fixed_params['profile_file']
    Bref         = fixed_params['Bref']
    Ma           = fixed_params['Ma']
    CDp          = fixed_params['CDp']

    AVL_path     = 'Results.avl'

    Yle = CalculateYle(Bref/2, sections)
    Xle = CalculateXle(chords, Yle)
    Zle = CalculateZle(dihed, Yle)

    Sref = CalculateSref(Yle, chords)
    Cref = CalculateMAC(Sref, Yle, chords)

    WriteAVLFile(
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

def main():
    config_path = "config.json"
    
    (
        file_name, wing_name, profile_name, 
        Starget, Bref, Ma, CLtarget, CDp, sections,
        x0, bounds
    ) = ReadConfigJSON(config_path)

    fixed_params = {
        'file_name': file_name,
        'wing_name': wing_name,
        'profile_file': profile_name,
        'Bref': Bref,
        'Ma': Ma,
        'CLtarget': CLtarget,
        'Starget':Starget,
        'CDp': CDp,
        'sections': sections
    }

    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')
    
    def AreaConstraint(params):
        chords   = params[1 : 1 + sections]
        Yle      = CalculateYle(Bref/2, sections)
        Sref     = CalculateSref(Yle, chords)
        return Sref - Starget
    

    # Nonlinear equality constraint: area ≈ Starget

    sref_constr = NonlinearConstraint(
        fun=AreaConstraint,
        lb=-0.01, # At 0.05 surface deviates too much:
        ub=0.01,  # Got values about 0.08-0.1 away from target
        jac="2-point"
    )
 
    """
    res = minimize(
        CostFunction, x0, args=(fixed_params),
        method='SLSQP', bounds=bounds, constraints=[constraints],
        options={"gtol": 1e-16, "disp": True, "eps": 1e-3}
        )"""
    
    res = differential_evolution(
        CostFunction,
        bounds=bounds,
        args=(fixed_params,),
        strategy="best1bin",
        popsize=20,
        tol=0.01,
        maxiter=20000,
        polish=True,
        constraints=(sref_constr,),
        workers=1,
        x0 = x0
    )

    logging.debug(res)

    WriteFinalResults(res.x, fixed_params)


if __name__ == '__main__':
    main()
