# Fargo Version of the CPT group

FargoCPT is a two-dimensional hydrodynamics + Nbody + particles code that is used to simulate protoplanetary disks with embedded planets and dust particles.

This code supports Linux (it runs in GitHub Codespaces!) and MacOS and is written in C++ using MPI and OpenMP for parallelization.

## Quickstart

The easiest way to run a simulation using FargoCPT is to launch a **codespace** from this repository.
To do this
1) Click the green Code button, and launch a codespace.
2) Wait a while (some minutes) until the codespace is generated.
3) Wait bit longer until the code is compiled (this should happen automatically).
4) Navigate to the Jupyter notebook located at [examples/100_Quickstart.ipynb](https://github.com/rometsch/fargocpt/blob/master/examples/100_Quickstart.ipynb).
5) Go through the notebook step by step.

This notebook will guide you through
- how to build the code,
- how to run a simulation, 
- describe the output structure 
- how to plot some results.
After completing this notebook, you likely know enough to run and analyze your own experiments.

To learn more about how to configure the code, please visit the [documentation](https://rometsch.github.io/fargocpt/).
It includes all simulation parameters and more examples to browse.
Alternatively, explore the `examples` directory.

Running the same notebooks locally, make sure to install the python dependencies (`pip install -r requirements.txt`) and go step by step through the notebooks.

## Running locally

To run locally on your computer, clone the last commit of the repository
``` bash
git clone https://github.com/rometsch/fargocpt --depth 1
```

## Description

This is the version of the [FARGO](http://fargo.in2p3.fr/-Legacy-archive-) code used by members of the Computational Physics Tübingen (CPT) group.

FargoCPT includes a number of improvements over the original FARGO code and is parallelized using a hybrid MPI-OpenMP scheme.

It is now used and maintained by Lucas Jordan and Thomas Rometsch with contributions by Tobias Moldenhauer (symmetric self-gravity) and Dennis Wehner (variable adiabatic index).

Credit for earlier version we built uopn go to:
- the original [FARGO](http://fargo.in2p3.fr/-Legacy-archive-) code by Frederic Masset
- [FARGO-ADSG](http://fargo.in2p3.fr/-FARGO-ADSG-) by [Clément Baruteau](http://clement.baruteau.free.fr/work/) who added a solver for the energy equation (AD=adiabatic) and self-gravity (SG)
- [Tobias Müller](https://twam.info) for adopting FARGO-ADSG to C++, [Giovanni Picogna](https://www.usm.uni-muenchen.de/people/picogna/index.html) for adding Langrangian particles, and other students in the CPT group.

## In development

This code is presented as is. And it will change over time.
If something does not work, please __open an issue on GitHub__.
We appreciate any feedback!

For changes, please see [ChangeLog.md](https://github.com/rometsch/fargocpt/blob/devel/ChangeLog.md)

## Usage

Start a simulation with one of the following commands:
```
./run_fargo -np NPROCS -nt NTHREADS {start/restart N/auto} setup/config.yml
```
where `NPROCS` is the number of MPI processes you want to launch and `NTHREADS` is the number of OpenMP threads per MPI process.

Using the `start` command begins a new simulation. 
Should the output directory already exists, the previous directory is backed up (`_bak` is appended to its name).

Using the `restart` command resumes a simulation that already exists in the output directory that is specified in the config file at output number `N`.

In `auto` mode, the simulation is restarted at the last available output if there are already some outputs, otherwise, a fresh simulation is started.

When omitting the `-np` and `-nt` options, the starting script tries to automatically use all available resources and determine the appropritate number of processes and threads to use for the given computer - one MPI process per numa node with as many OpenMP threads as there are cores per numa node. Hyperthreads are ignored.

`./run_fargo` is a script to call part of python module in `python-module` for calling the binary executable using `mpirun` with a couple of runtime options, e.g. for binding the processes to numa nodes.
We tested this setup for OpenMPI with versions 3 and 4.
For any other MPI implementation, please test the available options for hybrid parallelization using `mpirun` manually using
```bash
mpirun <your choice of options> ./bin/fargocpt_exe {start/restart N/auto} setup/config.yml
```

If you have the fargocpt python module installed (`pip install -r requirements`), you can also use the cli tool that is automatically installed with the python package.
Using this, you can run
```bash
fargocpt run --exe bin/fargocpt_exe -np NPROCS -nt NTHREADS {start/restart N/auto} setupfile.yml
```
Essentially, the `run_fargo` script is swapped out for `fargocpt run --exe ...`.
Please note the `--exe` option. You need to specify the path to your `fargocpt_exe` executable of you choice.
The program will complain, if it can't find it.

You can specify the executable path in three ways, where the first options take priority.
1. Using the `--exe` command line option (equivalently specifying the `exe=` keyword argument in python).
2. Specifying the path in `~/.config/fargocpt/config.yml` as key `exe_path: /path/to/fargocpt_exe`.
3. By placing `fargocpt_exe` into your `PATH` variable.

## Building the code

You likely want to install the python dependencies and the `fargocpt` python module.
Before or after making the code, run `pip install -r requirements`.

To build the code, navigate to the repository home and run

```bash
make -C src -j 4
```

This will compile the code in parallel using 4 processes. Increase this number at the risk of running out of memory.

The building process is managed by the makefile `src/makefile`.
Compile time flags are set in `src/makefile.defs` and the environment is selected and specified in `src/arch.defs`.

To get more information about the build process, run `make info`.

## Dependencies

Building and running FargoCPT requires the following dependencies:

- gcc / clang
- make
- git
- python3
- openmpi
- openmp (libomp)
- fftw (including lfftw3_mpi and lfftw3_omp)
- gsl

For anaconda users, please create a new envrioment and activate it.
Otherwise the mpi compiler installed with some packages might interfere with the compilation.

### Architecture definitions

FargoCPT depends on the three following libraries to be installed on your system:
- MPI
- OPENMP
- FFTW 
- GSL 

Their location can be specified using environment variables
- FARCOCPT_CC: select the c compiler (default mpicc)
- FARCOCPT_CXX: select the c++ compiler (default mpic++)
- FARCOCPT_CFLAGS: select additional CFLAGS (default empty)
- MPI_HOME: prefix of the MPI installation (default /usr)
- FFTW_HOME: prefix of the FFTW installation (default /usr)
- GSL_HOME: prefix of the GSL installation (default /usr)
- OMP_HOME: prefix of the GSL installation (default /usr)

If you installed MPI, FFTW, GSL, OMP through a package manager, chances are good that the defaults will work out of the box (e.g. for the Ubuntu example above).
The same can be expected for clusters that use the `module` framework. Then the `_HOME` variables should be set by loading the modules.

The `_HOME` variables should be the path that points to the directory that includes both the `include` and `lib` that contain the header files and shared libraries. To find them, search for, e.g., `mpi.h` and `libmpi.so` (use the `locate` or `find` commands). On MacOS, run `brew info {package name}`.


### Ubuntu/Debian

On Ubuntu (e.g. in a virtual maschine), run the following commands to get started.

``` bash
sudo apt-get install -y build-essential make
sudo apt-get install -y git
sudo apt-get install -y libopenmpi-dev
sudo apt-get install -y libgsl-dev
sudo apt-get install -y libfftw3-mpi-dev libfftw3-dev
sudo apt-get install -y python3
```

or as a one-liner
```bash
sudo apt-get install -y build-essential make git libopenmpi-dev  libgsl-dev libfftw3-mpi-dev libfftw3-dev python3
```

This has been tested on Ubuntu 20.04 and 22.04.

### MacOS

The following assumes you have a working `homebrew` installation. Tested on MacOS Ventura.

```zsh
brew install llvm # install an up to date compiler that supports openmp
brew install openmpi
brew install fftw
brew install gsl
brew install libomp
```


### Python

The wrapper for calling FargoCPT in parallel mode is written in Python3.
In case you want to run simulations in an automated fashion, you can import the `run_fargo` function from the `fargocpt.py` file as follows:

```python
import sys
sys.path.append("/path/to/the/repo/bin")
from fargocpt import run_fargo

N_procs = 2
N_OMP_threads = 8
run_fargo(N_procs, N_OMP_threads, ["start", "testconfig.yml"])
```

See `test/scaling` for an example use case.

## Docker

There is a `./docker/Dockerfile`, along with two `bash` scripts to build and run a docker image.

Call `./docker/build.sh` to build a docker image based on a Ubuntu 22.04 base.
Run the code with the `./docker/run.sh {mode} {setupfile} {outputdir}`.

Please note that OpenMPI can't bind memory to numa nodes when multiple processes are started. This effectively limits the docker image to be run on one numa node.
The `./docker/run.sh` script already takes this into account by launching only one process. It will be executed with as many threads as there are cores on one numa node.
This should be sufficient for many local applications.
On a cluster, you'll likely want to compile the code yourself.

The docker image assumes that the input files is located at `/simulation/setup.yml` inside the container and output files are written relative to `/simulation`.
Inside the container, the program can be called with the `fargocpt` command, which is a symlink to the python wrapper.


## Tests

There is a test suite to check that the physics modules work correctly.
Run it by executing `test/run_tests.sh`, it should run in a matter of a few minutes.
Each test outputs whether the test passed or failed.
Most of the test have a `readme.md` file with some description and produce a `plot.jpg` file to inspect the results.

- Shocktube (ideal and pvte equation of state)
- Steady state accretion disk (gas mass flow)
- Viscous spreading ring
- Viscous heating-cooling equilibrium temperature
- Irraditation-cooling equlibrium temperature
- Cold disk with ideal equation of state with and without planet
- Nbody integration: Kepler orbits
- Dust drift
- Dust diffusion
- Type I migration planet torque
- Flux-Limited Diffusion 1D
- Flux-Limited Diffusion 2D - direct test of the diffusion equation solver
- Self-gravity solver - comparison agains direct summation
- Planet orbiting a disk
- Steady-state accretion disk

## Open Source code used in this program

+ [REBOUND](https://github.com/hannorein/rebound) for N-body (GPL3)
+ [yaml-cpp](https://github.com/jbeder/yaml-cpp) for the config files (MIT License)
+ [units by LLNL](https://github.com/LLNL/units) for the handling of physical units (BSD 3-Clause)
+ [Opensimplex noise](https://github.com/deerel/OpenSimplexNoise), for randomization of the initial surface density profile ([The Unlicense
](https://github.com/deerel/OpenSimplexNoise/blob/master/UNLICENSE))
+ [jsf RNG](https://www.pcg-random.org/posts/bob-jenkins-small-prng-passes-practrand.html) for random number in dust diffusion and particle initialization. ([Source](https://gist.github.com/imneme/85cff47d4bad8de6bdeb671f9c76c814)) (MIT License)

