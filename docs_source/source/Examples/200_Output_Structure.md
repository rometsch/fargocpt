# FargoCPT Output Structure

This notebook introduces the structure of the simulation output directory and introduces a tool to investigate it.

We will use the data from the simulation in the quickstart example, so make sure you ran this beforehand.

## The output directory structure

The following cell will stop the notebook if the simulation has not been run yet.


```python
example_name = "100_quickstart"
example_dir = f"example_dirs/{example_name}"
import os
if not os.path.basename(os.getcwd()) == example_name:
    if os.path.exists(example_dir):
        os.chdir(example_dir)
        if not os.path.exists("output/out/snapshots/list.txt"):
            raise FileNotFoundError(f"No snapshots found. Please run the simulation inside the 100_Quickstart.ipynb notebook first!")
    else:
        raise FileNotFoundError(f"Please go through the 100_Quickstart.ipynb notebook first!")

repo_root = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
print(f"Current working directory: {os.getcwd()}")
print(f"Repository root directory: {repo_root}")
```

    Current working directory: /home/rometsch/repo/fargocpt/examples/example_dirs/100_quickstart
    Repository root directory: /home/rometsch/repo/fargocpt


Now, let's have a look at how the output directory is structured.


```python
!ls output/out
```

    constants.yml	      info1D.yml  monitor     units.yml
    dimensions.dat	      info2D.yml  parameters  used_rad.dat
    fargocpt_output_v1_4  logs	  snapshots


It contains some files describing general properties of the simulation, like 
- the version of the code (`fargocpt_output_v1_4`),
- the dimensions of the grid (`dimensions.dat` and `used_rad.dat`),
- the code units used (`units.dat`),
- and information on 1D output files (`*1D.info`) used to load the corresponding binary files.

There is also a `parameters` directory containing a copy of the setup file used for every start of the simulation.
This way you can easily track how often and at which snapshots you restarted a long run.
Here, we only have one copy.


```python
!ls output/out/parameters
```

    setup.yml


Next, there are snaphost directories, each containing a full snapshot of the system.
Each of these directories can be used to restart the simulation or start a new one.


```python
!ls output/out/snapshots
```

    0  1  10  2  3	4  5  6  7  8  9  list.txt  reference  timeSnapshot.dat


The `damping` directory contains a copy of the initial data which is used in the code e.g. for damping to the initial density inside of damping zones close to the boundaries. Copy this aswell, if you want to restart a simulation from a snapshot.

There is `timeSnapshot.dat` file which is a tab separated data file containing the time of the snapshot and the `list.txt` file which is simply a text file which has the number of each snapshot in a separate line. This is useful if you interact with the code using the command line. E.g.


```python
!tail -n 1 output/out/snapshots/list.txt
```

    10


is a clean way to get the number of the last snapshot.

The `snapshots` directory contains 
- the state variables of the hydro simulation (density, energy and velocities), 
- 1D output files, 
- a binary file for each planet and the `rebound.bin` for the state of the integrator (this is used for binary exact restarting)
- the `misc.bin` file which contains the state of the simulation system, e.g. the orientation of the coordinate system w.r.t. to an inertial frame and the last used CFL limited timestep,
- and a copy of the setup at the time of this snapshot.


```python
!ls output/out/snapshots/0
```

    config.yml    misc.bin	  rebound.bin  vazi1D.dat  vrad.dat
    energy1D.dat  nbody0.bin  Sigma1D.dat  vazi.dat
    energy.dat    nbody1.bin  Sigma.dat    vrad1D.dat


Finally, there is the `monitor` directory which contains monitor variables. These scalar variables are computed from the system state during the simulation and are written more often than the full snapshots.


```python
!ls output/out/monitor
```

    nbody0.dat  nbody1.dat	Quantities.dat	timeMonitor.dat  timestepLogging.dat


All files have a header that describes the colums and the units of the variables. The header can be automatically parsed.


```python
!head output/out/monitor/timeMonitor.dat
```

    # Time log for course output.
    #version: 0.1
    #variable: 0 | snapshot number | 1
    #variable: 1 | monitor number | 1
    #variable: 2 | time | 5.0225669513368811e+06 s
    # One monitor_timestep is 0.314000000000000001 (code) and 1577086.02271978068 (cgs).
    # Syntax: snapshot number <tab> monitor number <tab> time (cgs)
    0	0	0.0000000000000000e+00
    0	1	3.1400000000000000e-01
    0	2	6.2800000000000000e-01


- Each planet has its own file,
- disk quantities, e.g. the total mass, are stored in `Quantities.dat`,
- the output times of these fine grained monitor variables are stored in `timeMonitor.dat`, along with the corresponding snapshot number,
- and information about the CFL timestep and the ellapsed walltime can be found in `timestepLogging.dat`.

## Loading data

Let's inspect the monitor quantities that Fargo outputs.
Those are stored in the `monitor` directory within the output dir.

We'll use the `inspect_tab_file.py` tool, which helps navigating the tab separated output files.

Calling this tool with the `monitor/Quantities.dat` file, an overview of the available data is shown.


```python
!python3 $repo_root/Tools/inspect_tab_file.py output/out/monitor/Quantities.dat
```

    Available variables:
     0   snapshot number
     1   monitor number
     2   time
     3   mass
     4   radius
     5   angular momentum
     6   total energy
     7   internal energy
     8   kinematic energy
     9   potential energy
    10   radial kinetic energy
    11   azimuthal kinetic energy
    12   eccentricity
    13   periastron
    14   viscous dissipation
    15   luminosity
    16   pdivv
    17   inner boundary mass inflow
    18   inner boundary mass outflow
    19   outer boundary mass inflow
    20   outer boundary mass outflow
    21   wave damping inner mass creation
    22   wave damping inner mass removal
    23   wave damping outer mass creation
    24   wave damping outer mass removal
    25   density floor mass creation
    26   aspect ratio
    27   indirect term nbody x
    28   indirect term nbody y
    29   indirect term disk x
    30   indirect term disk y
    31   frame angle
    32   advection torque
    33   viscous torque
    34   gravitational torque



```python
!python3 $repo_root/Tools/inspect_tab_file.py output/out/monitor/Quantities.dat 2 3 --units kyr solMass | head
```

             0 kyr	  0.000349 solMass
         5e-05 kyr	  0.000349 solMass
      9.99e-05 kyr	  0.000349 solMass
       0.00015 kyr	  0.000349 solMass
        0.0002 kyr	  0.000349 solMass
       0.00025 kyr	  0.000349 solMass
        0.0003 kyr	  0.000349 solMass
       0.00035 kyr	  0.000349 solMass
        0.0004 kyr	  0.000349 solMass
       0.00045 kyr	  0.000349 solMass

