#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
OMP_NUM_THREADS=1 ../../fargo start dust_drift.yml 1> out.log 2>err.log
./calc_deviation.py