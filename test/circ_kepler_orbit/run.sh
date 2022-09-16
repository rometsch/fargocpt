#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
OMP_NUM_THREADS=1 ../../fargo start nbody_test.yml
./calc_deviation.py