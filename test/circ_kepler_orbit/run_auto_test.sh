#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
OMP_NUM_THREADS=1 ../../fargo start nbody_test.yml 1> out.log 2>err.log
./calc_deviation.py