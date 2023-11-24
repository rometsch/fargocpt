# FargoCPT Bootstrap Notebook

## What is this?
This notebook can be used to setup a FargoCPT simulation in an empty directory.

## Why is this useful?
This is useful when you want to run a simulation on a cluster which you can access using Jupyter notebooks and want to store a full copy of the code alongside your simulation outputs.

## Contents

It will do the following things:
- clone the repository from github
- compile the code
- give you a command to run the simulation
- print the available data
- plot the surface density

First we create a new directory and change to it.


```python
example_name = "300_bootstrap"
example_dir = f"example_dirs/{example_name}"
import os
if not os.path.basename(os.getcwd()) == example_name:
    !mkdir -p $example_dir
    os.chdir(example_dir)
repo_root = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
print(f"Current working directory: {os.getcwd()}")
print(f"Repository root directory: {repo_root}")
```

    Current working directory: /home/rometsch/repo/fargocpt/examples/example_dirs/300_bootstrap
    Repository root directory: /home/rometsch/repo/fargocpt


## Downloading the code

We will clone only the last commit of the code which is enough to run the simulation and faster to download.


```python
!git clone --depth 1 https://github.com/rometsch/fargocpt code
```

    Cloning into 'code'...
    remote: Enumerating objects: 749, done.[K
    remote: Counting objects: 100% (749/749), done.[K
    remote: Compressing objects: 100% (625/625), done.[K
    remote: Total 749 (delta 155), reused 459 (delta 88), pack-reused 0[K
    Receiving objects: 100% (749/749), 4.95 MiB | 4.37 MiB/s, done.
    Resolving deltas: 100% (155/155), done.


## Building the code

Make sure the code is built by running make again.


```python
%%timeit -n1 -r1
from sys import platform
if platform in ["linux", "darwin"]:
    !make -j 4 -C code/src > make.log
else:
    raise RuntimeError(f"Seems like you are not running MacOS or Linux but {platform}. This is unsupported. You are on your own, good luck!")
```

    26.7 s ± 0 ns per loop (mean ± std. dev. of 1 run, 1 loop each)


## Preparing a setup file

We'll take the example setup file from the examples directory and modify it in python.
If you want to create setup files for a parameter study, just copy the code and make your own setup creator script.


```python
configfile = "setup.yml"
!cp code/examples/config.yml $configfile
```

We'll use the `ruamel.yaml` package to read and write the setup file. This can be set up to preserve comments which is very useful if you want to trace your decisions later on.


```python
try:
    import ruamel.yaml
except ImportError:
    raise ImportError("Please install ruamel.yaml with `python3 -m pip install ruamel.yaml`")
yaml = ruamel.yaml.YAML()
with open(configfile, "r") as infile:
    config = yaml.load(infile)
```


```python
config["nbody"][1]["accretion efficiency"] = "2"
config["MonitorTimestep"] = 0.314 # monitor scalar files around every half orbit
config["Nmonitor"] = 20 # write a snapshot every orbit
config["Nsnapshots"] = 30 # wirte 10 snapshots
# use very low resolution by setting it to 2 cell per scaleheight, cps
del config["Nrad"]
del config["Naz"]
config["cps"] = 2
```


```python
with open(configfile, "w") as outfile:
    yaml.dump(config, outfile)
```

## Running the code

We can start fargo using the python interface, but this runs slower when started from within a Jupyter Notebook compared to being executed from without.
Even calling a python script using the shell magic "!" does not speed it up.
Calling a python script that does the same job from the command line does not have this issue.

If anyone knows why this is, please let me know!

For a production run, even for further testing, please open a terminal and run the output of the following cell:


```python
cwd = os.getcwd()
cmd = f"cd {cwd} && code/run_fargo -np 1 -nt 4 auto {configfile}"
print(cmd)
```

    cd /home/rometsch/repo/fargocpt/examples/example_dirs/300_bootstrap && code/run_fargo -np 1 -nt 4 auto setup.yml


If you plan to run on a cluster, put something similar into a run script to be queued in the queueing system.

Consider generating your queuing script here so you can just copy this notebook and setup a new simulation.

For the sake of this notebook, we just run the simluation here.

We'll use the python module from within the directory.


```python
import sys
sys.path = ["code/python"] + sys.path
from fargocpt import run
run(["start", configfile], np=2, nt=1, exe="code/bin/fargocpt_exe", detach=False)
```

    Running command: mpirun -np 2 --report-pid /tmp/tmpcqh1chwn -x OMP_NUM_THREADS=1 code/bin/fargocpt_exe start setup.yml
    fargo process pid 1403025
    
    [0] MPI rank #  0 runs as process 1403029
    [1] MPI rank #  1 runs as process 1403030
    [0] MPI rank #  0 OpenMP thread #  0 of  1 on cpt-kamino
    [1] MPI rank #  1 OpenMP thread #  0 of  1 on cpt-kamino
    [0] fargo: This file was compiled on Nov 14 2023, 15:52:04.
    [0] fargo: This version of FARGO used _GNU_SOURCE
    [0] fargo: This version of FARGO used NDEBUG. So no assertion checks!
    [0] Using parameter file setup.yml
    [0] Computing disk quantities within 5.00000e+00 L0 from coordinate center
    [0] BC: Inner composite = reflecting
    [0] BC: Outer composite = reflecting
    [0] BC: Sigma inner = zerogradient
    [0] BC: Sigma outer = zerogradient
    [0] BC: Energy inner = zerogradient
    [0] BC: Energy outer = zerogradient
    [0] BC: Vrad inner = reflecting
    [0] BC: Vrad outer = reflecting
    [0] BC: Vaz inner = keplerian
    [0] BC: Vaz outer = keplerian
    [0] DampingTimeFactor: 1.00000e-01 Outer damping time is computed at radius of 2.50000e+00
    [0] Damping VRadial to reference value at inner boundary.
    [0] Damping VRadial to reference value at outer boundary.
    [0] Damping VAzimuthal to reference value at inner boundary.
    [0] Damping VAzimuthal to reference value at outer boundary.
    [0] Damping SurfaceDensity to reference value at inner boundary.
    [0] Damping SurfaceDensity to reference value at outer boundary.
    [0] Damping Energy to reference value at inner boundary.
    [0] Damping Energy to reference value at outer boundary.
    [0] Radiative diffusion is disabled. Using fixed omega = 1.500000 with a maximum 50000 interations.
    [0] Indirect Term computed as effective Hydro center acceleratrion with shifting the Nbody system to the center.
    [0] Body force on gas computed via potential.
    [0] Using FARGO algorithm for azimuthal advection.
    [0] Using standard forward euler scheme for source terms.
    [0] Grid resolution set using cps = 2.000000
    [0] The grid has (Nrad, Naz) = (74, 251) cells with (1.994113, 1.997395) cps.
    [0] Computing scale height with respect to primary object.
    [0] Using isothermal equation of state. AdiabaticIndex = 1.400.
    [0] Viscosity is of alpha type with alpha = 1.000e-03
    [0] Defaulting to VanLeer flux limiter
    [0] Output information:
    [0]    Output directory: output/out/
    [0]     Number of files: 240
    [0]   Total output size: 0.00 GB
    [0]     Space Available: 31.18 GB
    [0] Initializing 1 RNGs per MPI process.
    [0] Warning : no `radii.dat' file found. Using default.
    [0] The first 1 planets are used to calculate the hydro frame center.
    [0] The mass of the planets used as hydro frame center is 1.000000e+00.
    [0] 2 planet(s) initialized.
    [0] Planet overview:
    [0] 
    [0]  #   | name                    | mass [m0]  | x [l0]     | y [l0]     | vx         | vy         |
    [0] -----+-------------------------+------------+------------+------------+------------+------------+
    [0]    0 | Star                    |          1 |          0 |         -0 |          0 |          0 |
    [0]    1 | Jupiter                 |  0.0009546033 |          1 |          0 |         -0 |   1.000477 |
    [0] 
    [0]  #   | e          | a          | T [t0]     | T [a]      | accreting  | Accretion Type |
    [0] -----+------------+------------+------------+------------+------------+----------------+
    [0]    0 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          0 |   No Accretion |
    [0]    1 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          2 |   Kley Accret. |
    [0] 
    [0]  #   | Temp [K]   | R [l0]     | irradiates | rampuptime |
    [0] -----+------------+------------+------------+------------+
    [0]    0 |       5778 |  0.0046505 |        yes |          0 |
    [0]    1 |          0 |  4.6505e-05 |         no |          0 |
    [0] 
    [0] Using Tscharnuter-Winkler (1979) artificial viscosity with C = 1.410000.
    [0] Artificial viscosity is used for dissipation.
    [0] Surface density factor: 2.50663
    [0] Tau factor: 0.5
    [0] Tau min: 0.01
    [0] Kappa factor: 1
    [0] Minimum temperature: 2.81162e-05 K = 3.00000e+00
    [0] Maximum temperature: 9.37206e+94 K = 1.00000e+100
    [0] Heating from viscous dissipation is enabled. Using a total factor of 1.
    [0] Cooling (beta) is disabled and reference temperature is floor. Using beta = 10.
    [0] Cooling (radiative) is enabled. Using a total factor of 1.
    [0] S-curve cooling is disabled. 
    [0] CFL parameter: 0.5
    [0] Opacity uses tables from Lin & Papaloizou, 1985
    [0] Particles are disabled.
    [0] Initializing Sigma(r) = 2.25093e-05 = 200 g cm^-2 * [r/(1 AU)]^(-0.5)
    [0] Total disk is mass is 0.000348841 = 6.9366e+29 g.
    [0] Writing output output/out/snapshots/0, Snapshot Number 0, Time 0.000000.
    [0] Writing output output/out/snapshots/reference, Snapshot Number 0, Time 0.000000.
    [0] Writing output output/out/snapshots/1, Snapshot Number 1, Time 6.280000.
    [0] Writing output output/out/snapshots/2, Snapshot Number 2, Time 12.560000.
    [0] Writing output output/out/snapshots/3, Snapshot Number 3, Time 18.840000.
    [0] Writing output output/out/snapshots/4, Snapshot Number 4, Time 25.120000.
    [0] Writing output output/out/snapshots/5, Snapshot Number 5, Time 31.400000.
    [0] Writing output output/out/snapshots/6, Snapshot Number 6, Time 37.680000.
    [0] Writing output output/out/snapshots/7, Snapshot Number 7, Time 43.960000.
    [0] Writing output output/out/snapshots/8, Snapshot Number 8, Time 50.240000.
    [0] Writing output output/out/snapshots/9, Snapshot Number 9, Time 56.520000.
    [0] Writing output output/out/snapshots/10, Snapshot Number 10, Time 62.800000.
    [0] Writing output output/out/snapshots/11, Snapshot Number 11, Time 69.080000.
    [0] Writing output output/out/snapshots/12, Snapshot Number 12, Time 75.360000.
    [0] Writing output output/out/snapshots/13, Snapshot Number 13, Time 81.640000.
    [0] Writing output output/out/snapshots/14, Snapshot Number 14, Time 87.920000.
    [0] Writing output output/out/snapshots/15, Snapshot Number 15, Time 94.200000.
    [0] Writing output output/out/snapshots/16, Snapshot Number 16, Time 100.480000.
    [0] Writing output output/out/snapshots/17, Snapshot Number 17, Time 106.760000.
    [0] Writing output output/out/snapshots/18, Snapshot Number 18, Time 113.040000.
    [0] Writing output output/out/snapshots/19, Snapshot Number 19, Time 119.320000.
    [0] Writing output output/out/snapshots/20, Snapshot Number 20, Time 125.600000.
    [0] Writing output output/out/snapshots/21, Snapshot Number 21, Time 131.880000.
    [0] Writing output output/out/snapshots/22, Snapshot Number 22, Time 138.160000.
    [0] Writing output output/out/snapshots/23, Snapshot Number 23, Time 144.440000.
    [0] Writing output output/out/snapshots/24, Snapshot Number 24, Time 150.720000.
    [0] Writing output output/out/snapshots/25, Snapshot Number 25, Time 157.000000.
    [0] Logging info: snapshot 25, monitor 514, hydrostep 4249, time inside simulation 161.619190, dt 4.524e-02, realtime 10.00 s, timeperstep 2.35 ms
    [0] Writing output output/out/snapshots/26, Snapshot Number 26, Time 163.280000.
    [0] Writing output output/out/snapshots/27, Snapshot Number 27, Time 169.560000.
    [0] Writing output output/out/snapshots/28, Snapshot Number 28, Time 175.840000.
    [0] Writing output output/out/snapshots/29, Snapshot Number 29, Time 182.120000.
    [0] Writing output output/out/snapshots/30, Snapshot Number 30, Time 188.400000.
    [0] -- Final: Total Hydrosteps 4909, Time 188.40, Walltime 11.56 seconds, Time per Step: 2.35 milliseconds





    0



## Data

Let's check whether there is data.


```python
from fargocpt import Loader
```


```python
l = Loader("output/out")
l
```




       Loader
    ====================
    | output_dir: output/out
    | snapshots: 0 ... 30
    | special_snapshots: ['reference']
    | snapshot_time: 0.0 5.02257e+06 s ... 188.4 5.02257e+06 s
    | monitor_number: 0 ... 600
    | units: Units
    | target_units = None
    | gas: Hydro
    | nbody: Nbody
    | params: Params
    | particles = None
    ====================




```python

```
