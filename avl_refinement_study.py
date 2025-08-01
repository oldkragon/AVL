import numpy as np
import subprocess
import logging
import re
import matplotlib.pyplot as plt
from multiprocessing import Pool


# Calculates Yle for each section with cosine distribution
def CalculateYle(b2:float, N:int):
    y_pan = b2*np.cos(np.pi/2 * np.linspace(0, N-1, N)/(N-1))
    return y_pan[::-1]


# Calculates chords for an elliptical wing
def CalculateChords(y_pan:np.ndarray, Sref:float, b2:float):
    c_root = 2 * Sref /(np.pi*b2)
    chords = c_root * np.sqrt(1-np.square(y_pan/b2))
    return chords


# Calculates the mean aerodynamic chord
def CalculateMAC(Sref:float, y_pan:np.ndarray, chords:list):
    if len(chords) != len(y_pan):
        raise ValueError("Chords and y_pan don't have the same length")
    chord_squared = np.square(np.array(chords))
    Mac = 2/Sref * np.trapz(chord_squared, y_pan)
    return Mac


# Calculates Xle for each point in y_pan with a sweep angle of dihed in DEGREES
def CalculateXle(chords:list, y_pan:np.ndarray, sweep=0.0):
    x_pan = np.array([chords[0] - c for c in chords])
    if sweep != 0.0:
        x_pan += y_pan* np.sin(np.pi/180 *sweep)
    return x_pan


# Calculates the wing surface
def CalculateSref(y_pan:np.ndarray, chords:list):
    if len(chords) != len(y_pan):
        raise ValueError("Chords and y_pan don't have the same length")
    Sref =  2 * np.trapz(chords, y_pan)
    return Sref


# Composes the .avl file for a surface
def WriteAVLFile(
        Ma:float, 
        sections:int,
        chords:list, 
        profile_file:str,
        AVL_path:str,
        Xle:np.ndarray, Yle:np.ndarray,
        Sref:float, Cref:float, Bref:float,
        num_pan_x:int, num_pan_y:int,
        num_pan_sect:int
        ):
    BaseAVLText =  f"""Refinement_study
{Ma:.2f}                   !   Mach
0     0     0.0       !   iYsym  iZsym  Zsym
{Sref:.4f} {Cref:.4f}  {Bref:.4f}       !   Sref   Cref   Bref   reference area, chord, span
0.0000   0.0000   0.0000       !   Xref   Yref   Zref   moment reference location (arb.)
0.012                 !   CDp
#
#==============================================================
#
SURFACE
Wing
{num_pan_x}  1.0  {num_pan_y}  1.0   ! Nchord   Cspace   Nspan  Sspace
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
    0.0000    0.0000     0.0000
"""
    for i in range(sections):
        SectionAVLText = f"""#
#--------------------------------------------------------------
#    Xle         Yle         Zle         chord       angle   Nspan  Sspace
SECTION
    {Xle[i]:.10f}         {Yle[i]:.10f}         0.0000         {chords[i]:.10f}        0.0000   {num_pan_sect}     -2

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
def RunAVL(Al_target:float, AVL_path:str, out_path:str):

    # Define the command sequence you want to send to AVL
    avl_commands = f"""
load {AVL_path}
oper
A
A {Al_target}
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
    

# Gets Cl, Cd and alpha from loads_file
def ParseResults(loads_file:str):
    with open(loads_file, 'r') as file:
        text = file.read()

    CL_match = re.search(r"CLtot\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)
    CD_match = re.search(r"CDtot\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)
    Alpha_match = re.search(r"Alpha\s*=\s*([-+]?[0-9]*\.?[0-9]+)", text)

    if not CL_match or not CD_match or not Alpha_match:
        raise ValueError("CL or CD not found in loads file.")
    
    CLtot = float(CL_match.group(1))
    CDtot = float(CD_match.group(1))
    Alpha = float(Alpha_match.group(1))

    return CLtot, CDtot, Alpha


def AVLCase(
        
        profile_file, 
        avl_path, loads_path, 
        Ma, Al_target,
        Starget, Bref, 
        Num_sect,
        num_pan_x, num_pan_sect):
    
    Yle    = CalculateYle(Bref/2, Num_sect)
    chords = CalculateChords(Yle, Starget, Bref/2)
    Sref   = CalculateSref(Yle, chords)
    Cref   = CalculateMAC(Sref, Yle, chords)
    Xle    = CalculateXle(chords, Yle)

    loads_filename = f"temp/loads_x_{num_pan_x}_y_{num_pan_sect}_sect_{Num_sect}.txt"
    avl_filename   = f"temp/input_x_{num_pan_x}_y_{num_pan_sect}_sect_{Num_sect}.avl"

    num_pan_y = Num_sect * num_pan_sect

    WriteAVLFile(
        Ma,
        Num_sect,
        chords,
        profile_file,
        avl_filename,
        Xle, Yle,
        Sref, Cref, Bref,
        num_pan_x, num_pan_y,
        num_pan_sect
    )

    RunAVL(
        Al_target,
        avl_filename,
        loads_filename
    )

    Cl, Cd, alpha = ParseResults(loads_filename)
    return Cl, Cd, alpha


def XpanRefinementStudyParallel(
        max_num_pan_x,
        profile_file = 's9037opt5.dat',
        avl_path     = 'temp/test.avl',
        loads_path   = 'temp/test.txt',
        Ma           = 0.053,
        Al_target    = 0.5,
        Starget      = 0.85,
        Bref         = 3,
        Num_sect     = 10,
        num_pan_y    = 40,
        num_pan_sect = 4
        ):
    
    num_pans_x = np.linspace(1, max_num_pan_x, max_num_pan_x)

    args_list = [(profile_file, avl_path, loads_path, Ma, Al_target,
                  Starget, Bref, Num_sect, num_pan_x, num_pan_sect)
                 for num_pan_x in num_pans_x]

    with Pool(processes=15) as pool:
        results = pool.starmap(AVLCase, args_list)

    Cl_list, Cd_list, Alpha_list = zip(*results)
    return Cl_list, Cd_list, Alpha_list


def YpanRefinementStudyParallel(
        max_num_pan_sect,
        profile_file = 's9037opt5.dat',
        avl_path     = 'temp/test.avl',
        loads_path   = 'temp/test.txt',
        Ma           = 0.053,
        Al_target    = 0.5,
        Starget      = 0.85,
        Bref         = 3,
        Num_sect     = 10,
        num_pan_x    = 10,
        ):
    
    num_pans_sect = np.linspace(2, max_num_pan_sect, max_num_pan_sect-1, dtype=int)

    args_list = [(profile_file, avl_path, loads_path, Ma, Al_target,
                  Starget, Bref, Num_sect, num_pan_x, num_pan_sect)
                 for num_pan_sect in num_pans_sect]

    with Pool(processes=15) as pool:
        results = pool.starmap(AVLCase, args_list)

    Cl_list, Cd_list, Alpha_list = zip(*results)
    return Cl_list, Cd_list, Alpha_list


def NumSectRefinementStudyParallel(
        max_num_sect,
        profile_file = 's9037opt5.dat',
        avl_path     = 'temp/test.avl',
        loads_path   = 'temp/test.txt',
        Ma           = 0.053,
        Al_target    = 0.5,
        Starget      = 0.85,
        Bref         = 3,
        num_pan_x    = 10,
        num_pan_sect = 4
        ):
    
    num_sects = np.linspace(2, max_num_sect, max_num_sect-1, dtype=int)

    args_list = [(profile_file, avl_path, loads_path, Ma, Al_target,
                  Starget, Bref, num_sect, num_pan_x, num_pan_sect)
                 for num_sect in num_sects]

    with Pool(processes=15) as pool:
        results = pool.starmap(AVLCase, args_list)

    Cl_list, Cd_list, Alpha_list = zip(*results)
    return Cl_list, Cd_list, Alpha_list


# Old single-threaded functions, use Parallel version instead
def XpanRefinementStudy(
        max_num_pan_x,
        profile_file = 's9037opt5.dat',
        avl_path     = 'temp/test.avl',
        loads_path   = 'temp/test.txt',
        Ma           = 0.053,
        Cl_target    = 0.415,
        Starget      = 0.85,
        Bref         = 3,
        Num_sect     = 10,
        num_pan_y    = 40,
        num_pan_sect = 4
        ):
    
    Cl_list    = []
    Cd_list    = []
    Alpha_list = []

    for num_pan_x in range(max_num_pan_x):
        Yle    = CalculateYle(Bref/2, Num_sect)
        chords = CalculateChords(Yle, Starget, Bref/2)
        Sref   = CalculateSref(Yle, chords)
        Cref   = CalculateMAC(Sref, Yle, chords)
        Xle    = CalculateXle(chords, Yle)

        WriteAVLFile(
            Ma,
            Num_sect,
            chords,
            profile_file,
            avl_path,
            Xle, Yle,
            Sref, Cref, Bref,
            num_pan_x, num_pan_y,
            num_pan_sect
        )

        RunAVL(
            Cl_target,
            avl_path,
            loads_path
        )

        (Cl, Cd, alpha) = ParseResults(loads_path)

        Cl_list.append(Cl)
        Cd_list.append(Cd)
        Alpha_list.append(alpha)
    
    return Cl_list, Cd_list, Alpha_list


def YpanRefinementStudy(
        max_num_pan_sect,
        profile_file = 's9037opt5.dat',
        avl_path     = 'temp/test.avl',
        loads_path   = 'temp/test.txt',
        Ma           = 0.053,
        Cl_target    = 0.415,
        Starget      = 0.85,
        Bref         = 3,
        Num_sect     = 10,
        num_pan_x    = 10,
        ):
    
    Cl_list    = []
    Cd_list    = []
    Alpha_list = []

    for num_pan_sect in range(max_num_pan_sect):
        num_pan_y = num_pan_sect * Num_sect

        Yle    = CalculateYle(Bref/2, Num_sect)
        chords = CalculateChords(Yle, Starget, Bref/2)
        Sref   = CalculateSref(Yle, chords)
        Cref   = CalculateMAC(Sref, Yle, chords)
        Xle    = CalculateXle(chords, Yle)

        WriteAVLFile(
            Ma,
            Num_sect,
            chords,
            profile_file,
            avl_path,
            Xle, Yle,
            Sref, Cref, Bref,
            num_pan_x, num_pan_y,
            num_pan_sect
        )

        RunAVL(
            Cl_target,
            avl_path,
            loads_path
        )

        (Cl, Cd, alpha) = ParseResults(loads_path)

        Cl_list.append(Cl)
        Cd_list.append(Cd)
        Alpha_list.append(alpha)
    
    return Cl_list, Cd_list, Alpha_list


def NumSectRefinementStudy(
        max_num_sect,
        profile_file = 's9037opt5.dat',
        avl_path     = 'temp/test.avl',
        loads_path   = 'temp/test.txt',
        Ma           = 0.053,
        Al_target    = 0.5,
        Starget      = 0.85,
        Bref         = 3,
        num_pan_x    = 10,
        num_pan_sect = 4
        ):
    
    Cl_list    = []
    Cd_list    = []
    Alpha_list = []

    for num_sect in range(2, max_num_sect):
        num_pan_y = num_sect * num_pan_sect

        Yle    = CalculateYle(Bref/2, num_sect)
        chords = CalculateChords(Yle, Starget, Bref/2)
        Sref = CalculateSref(Yle, chords)
        Cref   = CalculateMAC(Sref, Yle, chords)
        Xle    = CalculateXle(chords, Yle)

        WriteAVLFile(
            Ma,
            num_sect,
            chords,
            profile_file,
            avl_path,
            Xle, Yle,
            Sref, Cref, Bref,
            num_pan_x, num_pan_y,
            num_pan_sect
        )

        RunAVL(
            Al_target,
            avl_path,
            loads_path
        )

        (Cl, Cd, alpha) = ParseResults(loads_path)

        Cl_list.append(Cl)
        Cd_list.append(Cd)
        Alpha_list.append(alpha)
    
    return Cl_list, Cd_list, Alpha_list


def main():
    max_num_pan_x = 20
    x = np.linspace(1, max_num_pan_x, max_num_pan_x)
    max_num_pan_sect = 20
    y = np.linspace(2, max_num_pan_sect, max_num_pan_sect-1)
    max_num_sect = 20
    z = np.linspace(2, max_num_sect, max_num_sect-1)

    study_name = input('Insert study name')
 
    (Cl_x, Cd_x, Alpha_x) = XpanRefinementStudyParallel(max_num_pan_x)
    (Cl_y, Cd_y, Alpha_y) = YpanRefinementStudyParallel(max_num_pan_sect)
    (Cl_z, Cd_z, Alpha_z) = NumSectRefinementStudyParallel(max_num_sect)

    plt.figure(1)

    plt.subplot(1, 3, 1)
    plt.plot(x, Cl_x)
    plt.title('Cl - num_pan_x')
    plt.grid(True)

    plt.subplot(1, 3, 2)
    plt.plot(x, Cd_x)
    plt.title('Cd - num_pan_x')
    plt.grid(True)

    plt.subplot(1, 3, 3)
    plt.plot(x, Alpha_x)
    plt.title('Alpha - num_pan_x')
    plt.grid(True)

    plt.tight_layout()

    plt.figure(2)

    plt.subplot(1, 3, 1)
    plt.plot(y, Cl_y)
    plt.title('Cl - num_pan_y')
    plt.grid(True)

    plt.subplot(1, 3, 2)
    plt.plot(y, Cd_y)
    plt.title('Cd - num_pan_y')
    plt.grid(True)

    plt.subplot(1, 3, 3)
    plt.plot(y, Alpha_y)
    plt.title('Alpha - num_pan_y')
    plt.grid(True)

    plt.tight_layout()

    plt.figure(3)

    plt.subplot(1, 3, 1)
    plt.plot(z, Cl_z)
    plt.title('Cl - num_sect')
    plt.grid(True)

    plt.subplot(1, 3, 2)
    plt.plot(z, Cd_z)
    plt.title('Cd - num_sect')
    plt.grid(True)

    plt.subplot(1, 3, 3)
    plt.plot(z, Alpha_z)
    plt.title('Alpha - num_sect')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    with open('Study.txt', 'a') as file:
        file.write(
            f'''Results {study_name}:
    Panels X CDs:
{Cd_x}

    Panels Y CDs:
{Cd_y}

    Sections CDs:
{Cd_z}

'''
    )


if __name__ == '__main__':
    main()