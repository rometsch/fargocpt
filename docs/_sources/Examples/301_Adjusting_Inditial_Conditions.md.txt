# Adjusting initial conditions

This notebook teaches you how to change the initial conditions of the disk.

First we create a new directory and change to it.


```python
example_name = "301_Adjusting_Inditial_Conditions"
example_dir = f"example_dirs/{example_name}"
import os
if not os.path.basename(os.getcwd()) == example_name:
    !mkdir -p $example_dir
    os.chdir(example_dir)
repo_root = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
print(f"Current working directory: {os.getcwd()}")
print(f"Repository root directory: {repo_root}")
```

    Current working directory: /home/rometsch/repo/fargocpt/examples/example_dirs/301_Adjusting_Inditial_Conditions
    Repository root directory: /home/rometsch/repo/fargocpt


## Make sure the code is built by running make again.


```python
%%timeit -n1 -r1
from sys import platform
if platform in ["linux", "darwin"]:
    !make -j 4 -C $repo_root/src > make.log
else:
    raise RuntimeError(f"Seems like you are not running MacOS or Linux but {platform}. This is unsupported. You are on your own, good luck!")
```

    make: *** ../../src: No such file or directory.  Stop.
    110 ms ± 0 ns per loop (mean ± std. dev. of 1 run, 1 loop each)


## Preparing a setup file

We'll take the example setup file from the examples directory and modify it in python.
If you want to create setup files for a parameter study, just copy the code and make your own setup creator script.


```python
configfile = "setup.yml"
!cp $repo_root/examples/config.yml $configfile
```

We'll use the `ruamel.yaml` package to read and write the setup file. This can be set up to preserve comments which is very useful if you want to trace your decisions later on.


```python
try:
    import ruamel.yaml
except ImportError:
    raise ImportError("Please install ruamel.yaml with `python3 -m pip install ruamel.yaml`")
yamlparser = ruamel.yaml.YAML()
with open(configfile, "r") as infile:
    config = yamlparser.load(infile)
```


```python
# we don't need to run long to inspect the initial conditions
config["MonitorTimestep"] = 0.314 # monitor scalar files around every half orbit
config["Nmonitor"] = 1 # write a snapshot every orbit
config["Nsnapshots"] = 1 # wirte 100 snapshots

# use very low resolution by setting it to 2 cell per scaleheight, cps
config["cps"] = 2

# set initial conditions of the powerlaw disk
config["SigmaSlope"] = 1.5 # Sigma(r) = Sigma0 * r**(-SigmaSlope)
config["Sigma0"] = "100 g/cm2" # we can use units here

# also, lets set a different temperature profile via the aspect ratio
config["AspectRatio"] = 0.1

# write out temperature and aspect ratio, so we can plot it
config["WriteTemperature"] = True
config["WriteAspectRatio"] = True

with open(configfile, "w") as outfile:
    yamlparser.dump(config, outfile)
```

## Run the simulation


```python
from fargocpt import run
run(["start", configfile], np=2, nt=1, exe=repo_root+"/bin/fargocpt_exe", detach=False)
```

    Running command: mpirun -np 2 --report-pid /tmp/tmp34swwsrv -x OMP_NUM_THREADS=1 /home/rometsch/repo/fargocpt/bin/fargocpt_exe start setup.yml
    fargo process pid 1410280
    
    [0] MPI rank #  0 runs as process 1410284
    [1] MPI rank #  1 runs as process 1410285
    [1] MPI rank #  1 OpenMP thread #  0 of  1 on cpt-kamino
    [0] MPI rank #  0 OpenMP thread #  0 of  1 on cpt-kamino
    [0] fargo: This file was compiled on Nov 14 2023, 12:56:40.
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
    [0] Cps is set, overwriting Nrad and Naz!
    [0] Grid resolution set using cps = 2.000000
    [0] The grid has (Nrad, Naz) = (38, 127) cells with (2.023980, 2.021268) cps.
    [0] Computing scale height with respect to primary object.
    [0] Using isothermal equation of state. AdiabaticIndex = 1.400.
    [0] Viscosity is of alpha type with alpha = 1.000e-03
    [0] Defaulting to VanLeer flux limiter
    [0] Output information:
    [0]    Output directory: output/out/
    [0]     Number of files: 12
    [0]   Total output size: 0.00 GB
    [0]     Space Available: 31.24 GB
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
    [0]    1 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          0 |   No Accretion |
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
    [0] Initializing Sigma(r) = 1.12546e-05 = 100 g cm^-2 * [r/(1 AU)]^(-1.5)
    [0] Total disk is mass is 0.000134118 = 2.6669e+29 g.
    [0] Writing output output/out/snapshots/0, Snapshot Number 0, Time 0.000000.
    [0] Writing output output/out/snapshots/reference, Snapshot Number 0, Time 0.000000.
    [0] Writing output output/out/snapshots/1, Snapshot Number 1, Time 0.314000.
    [0] -- Final: Total Hydrosteps 5, Time 0.31, Walltime 0.02 seconds, Time per Step: 3.64 milliseconds





    0



Following is an overview widget for the simulation. You can use the slider to scrub through the different snapshots.

Let's see which variables we have.


```python
from fargocpt import Loader
l = Loader("output/out")
l.gas.vars2D
```




       Vars2D
    ====================
    | output_dir: output/out
    | target_units= None
    | grid: Grid
    | var_names:
    |   Sigma
    |   vrad
    |   vazi
    |   energy
    |   Temperature
    |   aspectratio
    ====================



Run the next cell again to refresh the snapshot list.


```python
# %matplotlib widget
from fargocpt import Overview
overview = Overview("output/out/", 
                    vars=["1:Sigma:minmax",
                          "1:Temperature"])
overview.create();
```


    
![png](301_Adjusting_Inditial_Conditions_files/301_Adjusting_Inditial_Conditions_17_0.png)
    



```python

```
