# AVL scripts
This folder contains python scripts that call Athena Vortex Lattice for 3D wing optimization.

## AVL_script.py
This is the main script, that performs the wing optimization using *AVL* and scipy's *differential_evolution*.
To work correctly AVL needs to be added to the PATH environment variable.

To add AVL to PATH copy the path to the folder that contains avl.exe, then go into settings -> edit the system environment variables -> advanced -> Environment Variables -> Path -> edit -> new and paste the path.

To check that avl has been correctly added to path open a terminal and type avl, you should see something like " Athena Vortex Lattice  Program      Version  X.XX" 

### Input files

The script takes two file to run correctly:
- First a .dat file containing the coordinates of the airfoil used, formatted in XFOIL's style (the first row is the airfoil's name)
- Second a .json file containing the parameters for the optimization, as of right now this file MUST be called "config.json" and be placed in the same directory as avl_script.py

### Config file options
- file_name : this will be the first row of every file fed to AVL, to keep track of which optimization is being done.
- wing_name : in the file fed to AVL is the name given to the surface.
- profile_file : specify the relative or absolute path to the .dat file for the airfoil
- num_sect : the number of sections in which the wing will be divided by AVL
- CDp : default profile drag coefficient added to geometry
- opt_points :  a list of the optimization points, each point takes the following parameters
	-  name : not read by the script, used as a comment to distinguish points
	- op_mode : either 'spec_cl' where you provide S Cl or 'spec_al' where you provide alpha
	- op_point : the value to provide according to the "op_mode"
	- Ma : Mach's number
	- weight : you can assign different weights to each point, it is recommended, but not necessary, to have them add up to one
- span, twist and dihed : fields to provide bounds for the optimization, each has to include a sub-field "initial" for an initial guess and "bounds" to limit the region to explore
- sections : a list of bounds and initial guesses for the chords at each section, with the following parameters
	- name : not read by the script, used as a comment to distinguish points
	- chord:
		- bounds
		- initial
As an example:
```json
{
  "file_name": "Olivia",
  "wing_name": "Wing",
  "profile_file": "s9037opt5.dat",
  "num_sect": 10,
  "Starget": 0.85,
  "CDp": 0.012,
  "opt_points": [
    {
      "name": "with payload",
      "op_mode": "spec_cl",
      "op_point": 0.415,
      "Ma": 0.053,
      "u_inf": 18.0,
      "weight": 0.7
    }
  ],
  "span": {
    "bounds": [1.0, 5.0],
    "initial": 3.0
  },
  "dihed": {
    "bounds": [0.0, 2.0],
    "initial": 0.1
  },
  "twist": {
    "bounds": [-2.0, 0.5],
    "initial": -0.1
  },
  "sections": [
    {
      "name": 1,
      "chord": {
        "bounds": [0.05, 0.5],
        "initial": 0.284
      }
    }
  ]
}
``` 