# Fargo Version of the CPT group

## Description

This is the version of the [FARGO](http://fargo.in2p3.fr/-Legacy-archive-) code used by some members of the Computational Physics Tübingen (CPT) group.

It is now used and maintained by Lucas Jordan and Thomas Rometsch with contributions by Tobias Moldenhauer (symmetric self-gravity) and Dennis Wehner (variable adiabatic index).

Credit for earlier version we built uopn go to:
- the original [FARGO](http://fargo.in2p3.fr/-Legacy-archive-) code by Frederic Masset
- [FARGO-ADSG](http://fargo.in2p3.fr/-FARGO-ADSG-) by [Clément Baruteau](http://clement.baruteau.free.fr/work/) who added a solver for the energy equation (AD=adiabatic) and self-gravity (SG)
- [Tobias Müller](https://twam.info) for adopting FARGO-ADSG to C++, [Giovanni Picogna](https://www.usm.uni-muenchen.de/people/picogna/index.html) for adding Langrangian particles, and other students in the CPT group.

## Warning!

This version is still in development.
Most physics modules are successfully tested, however some modules (dust diffusion, self-gravity) are not yet fully tested.

The tested features of the code include:
- hydrodynamics part including the FARGO algorithm (shocktube, viscous spreading ring)
- Nbody dynamics and interaction with gas (planet torque, barycenter test)
- irradiation
- dust drift

## In development

This code is presented as is. And it will change over time.
The coming weeks (last weeks of 2022 and first weeks of 2023), we plan to add more tests, add documentation and examples.
Changes to the input and output format are also possible, though the goal for the next month is to have a version 1.0.

## Usage

Start a simulation with one of the following commands:
```
mpirun -n NPROC ./fargo start setup/config.yml
mpirun -n NPROC ./fargo restart 5 setup/config.yml
mpirun -n NPROC ./fargo auto setup/config.yml
```
where `NPROC` is the number of processors you want to use.

The first line with the `start` command begins a new simulation. 
Should the output directory already exists, the previous directory is backed up (`_bak` is appended to its name).

The second line is used to `restart` a simulation that already exists in the output directory that is specified in the config file at output number 5.

The third line starts fargo in `auto` mode. In this mode, the simulation is restarted at the last available output if there are already some outputs, otherwise, a fresh simulation is started.


## Building the code

The building process is managed by the makefile `src/makefile`.
Compile time flags are set in `src/makefile.defs` and the environment is selected and specified in `src/arch.defs`.

To get more information about the build process, run `make -m`

### Architecture definitions

FargoCPT depends on the three following libraries to be installed on your system:
- MPI
- FFTW 
- GSL 

Their location can be specified using environment variables
- FARCOCPT_CC: select the c compiler (default mpicc)
- FARCOCPT_CXX: select the c++ compiler (default mpic++)
- FARCOCPT_CFLAGS: select additional CFLAGS (default empty)
- MPI_HOME: prefix of the MPI installation (default /usr)
- FFTW_HOME: prefix of the FFTW installation (default /usr)
- GSL_HOME: prefix of the GSL installation (default /usr)

If you installed MPI, FFTW, and GSL through a package manager, chances are good that the defaults will work out of the box.
The same can be expected for clusters that use the `module` framework. Then the `_HOME` variables should be set by loading the modules.

The `_HOME` variables should be the path that points to the directory that includes both the `include` and `lib` that contain the header files and shared libraries. To find them, search for, e.g., `mpi.h` and `libmpi.so` (use the `locate` or `find` commands).

## TODOs

+ move subkeplerian boundary condition call from artificalviscosity.cpp to boundary conditions / simulation.cpp 

## Tests

+ Hydro
  + Shocktube
+ Viscosity
  + Viscous spreading ring
+ Heating/Cooling
  + passive disk -> power law profile
  + beta cooling background
  + beta cooling aspect ratio
+ Dust
  + Dust drift
  + Dust diffusion
+ N-body
  + migrating low mass planet
+ self-gravity
TBD
+ program
  + valgrind run
    + ignore OpenMPI errors
      + ignore rebound padding (Syscall param write(buf) points to uninitialized byte(s))


## Open Source code used in this program

+ [REBOUND](https://github.com/hannorein/rebound) for N-body (GPL3)
+ [yaml-cpp](https://github.com/jbeder/yaml-cpp) for the config files (MIT License)
+ [units by LLNL](https://github.com/LLNL/units) for the handling of physical units (BSD 3-Clause)