#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
../../run_fargo -nt 1 start nbody_test.yml 1> sim.log 2> sim.err
./plot_nbody.py
./calc_deviation.py