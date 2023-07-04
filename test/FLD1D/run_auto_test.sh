#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR

# we need to modify the build to enable a constant flux limiter
# ./modify_build.sh > modify_build.log


sed -i 's/RadiativeDiffusionVariable:.*$/RadiativeDiffusionVariable: temperature/' setup.yml 1> out_energy.log 2>err_energy.log
../../run_fargo -nt 2 -np 1 start setup.yml 1> out_temperature.log 2>err_temperature.log
./calc_deviation.py temperature
./plot_overview.py output/out overview_temperature.jpg

sed -i 's/RadiativeDiffusionVariable:.*$/RadiativeDiffusionVariable: energy/' setup.yml 
../../run_fargo -nt 2 -np 1 start setup.yml 1> out_energy.log 2>err_energy.log
./calc_deviation.py energy
./plot_overview.py output/out overview_energy.jpg



# revert the changes to the code
# ./revert_build.sh > revert_build.log
