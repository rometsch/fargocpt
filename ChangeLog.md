# Changelog of FargoCPT

Changes are in inverse chronological order, latest changes come first.

Prior to version 1.3.0, there was no versioning, so we used the latest output version (1.3) and made it the latest code version.

Syntax of the versions are: major, minor, revision
- major: major restructuring of code or config files or output
- minor: new features, new options, new config file options, renamed parameters, new or changed single output files
- revision: bug fixes

## Version 1.4
- rename 'vtheta' and 'vphi' both to 'vazi'. The output file of azimuthal velocity is now called 'vazi.dat' instead of 'vtheta.dat'. You might have to rename files to restart data.
- executable file renamed from `fargocpt` to `fargocpt_exe` to allow python cli to be called `fargocpt`
- added installable python package (`pip install <repo path>/python`) which installs a command line interface called `fargocpt`. Start a simulation with `fargocpt run` and check the data with `fargocpt data <output dir>`
- `planet{n}.dat/bin` files renamed to `nbody{n-1}.dat/bin`, thus `nbody0.dat` instead of `planet1.dat`. This requires renaming files to restart older data.
- nbody indices start at 0 now instead of at 1.

## Version 1.3
- the code now has a version string, run ./bin/fargocpt --version
- there is now a `parameters.md` file with a table of all parameters.- changed the name of some parameters including Ntot->Nsnapshot, Ninterm->Nmonitor, DT->MonitorTimestep, run ./Tools/param_names/replace_parameter_names.py on a yml config file to update the names.
- code now exits in case of unknown parameters in the setup file. This help spotting typos.
- info about particle data is now in the `infoParticles.yml` file
- info about 1D output is now in the `info1D.yml` file
- copies of log written to stdout are written to a output and error file in the log directory in the output directory
- backtraces are written to a separate file for each process
- constants file in yaml containing all physical constants in the code
- units file in yaml containing all units in the code including code unit value in cgs and the cgs unit symbols
- Radiative Diffusion Tolerance parameter is now an absolute quantity, if only a number is supplied, its multiplied by the Temperature unit, otherwise its parsed.
- there is a file `info2D.yml` which describes the 2D output variables including their array size and units
- self gravity acceleration can now be written out to the snapshots using `WriteSGAccelRad` and `WriteSGAccelAzi` parameters
- Beta cooling includes a new `floor` mode in which the temperature is relaxed to the floor value. Can be configured using the `CoolingBetaReference` parameter which can be `zero` (to zero Kelvin), `reference` (to the reference/damping values), `model` (to the temperature from aspect ratio). This will break your old config.
- Cooling from the disk surface: the config parameter `CoolingRadiativeLocal` has been replaced by `SurfaceCooling` which can take the value `no|false|off`, `thermal` (the existing mode), or `scurve` (new mode).
- there is a file `info2D.yml` which describes the 2D output variables including their array size and units
- self gravity acceleration can now be written out to the snapshots using `WriteSGAccelRad` and `WriteSGAccelAzi` parameter
- bigplanet.dat files are now called planet.dat files, now fargocpt output v1.2
- restructured the output directory into snapshots and monitor directories. Each snapshot dir contains all information to restart the simulation.
- replaced the default executable to run, now run_fargo starts the simulations and calls bin/fargocpt and handles the parallel config
- simplify build process and remove FARGO_ARCH, use specific env variables now (see arch.defs)
- backtraces are written to file by default
- irradiation is now controlled by setting a temperature > 0 in the nbody config rather than the "HeatingStar" flag
- removed default star, central object needs to be declared as a planet now
- OUTPUTDIR is now in output namespace and called outdir
- physical constants (gravitational constant, kB, atomic mass unit, Planck constant) from LLNL units library which follows the NIST constants and the 2019 redefinition of SI
- stellar heating rampup time in units of T0 rather than DT
- beta cooling rampup time in units of T0
- thinknesssmoothing factor 0.6 by default (was 0.0)
- Changed tau_factor default value from 1 to 0.5.