#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
../../run_fargo -nt 1 start dust_drift.yml 1> out.log 2>err.log
python3 calc_deviation.py
python3 plot_drift.py
python3 plot_trajectory.py