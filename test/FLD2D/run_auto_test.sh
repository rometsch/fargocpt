#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR

# create a sceleton output first

../../run_fargo -nt 1 -N 0 start setup.yml 1> out.log 2>err.log
./create_input.py

# ../../run_fargo -nt 1 -v restart 0 setup.yml # 1> out.log 2>err.log
mpirun -x OMP_NUM_THREADS=1 -n 1 ../../bin/fargocpt -v restart 0 setup.yml