#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR

# we need to modify the build to enable a constant flux limiter
# ./modify_build.sh > modify_build.log

../../run_fargo -nt 2 -np 1 start setup.yml 1> out.log 2>err.log
./calc_deviation.py

# revert the changes to the code
# ./revert_build.sh > revert_build.log
