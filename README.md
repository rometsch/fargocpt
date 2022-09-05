# Fargo Version of the CPT group

## Description

This is the version of the [FARGO](http://fargo.in2p3.fr/-Legacy-archive-) code used by some members of the Computational Physics Tübingen group.
It was modified quite substantially by Tobias Müller and is now used and maintained by Lucas Jordan and Thomas Rometsch.


## Building the code

The building process is managed by the makefile `src/makefile`.
Compile time flags are set in `src/makefile.defs` and the environment in selected and specified in `src/arch.defs`.

### 


## Tests

+ Hydro
    + Shocktube
+ Viscosity
    + Viscous spreading ring
+ Heating/Cooling
    + passive disk -> power law profile
    + betacooling background
    + betacooling aspect ratio
+ Dust
    + Dust drift
    + Dust diffusion
+ Nbody
    + migrating low mass planet
+ self-gravity
    + tbd
+ program
    + valgrind run
        + ignore openmpi errors
        + ignore rebound padding (Syscall param write(buf) points to uninitialised byte(s))


### Open Source code used in this program

+ [REBOUND](https://github.com/hannorein/rebound) for NBody (GPL3)
+ [yaml-cpp](https://github.com/jbeder/yaml-cpp) for the config files (MIT License)
+ [units by LLNL](https://github.com/LLNL/units) for the handling of physical units (BSD 3-Clause)