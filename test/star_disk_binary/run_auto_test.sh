#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
OMP_NUM_THREADS=1 ../../fargo start leapfrog.yml 1> LF_disk_at_center.log 2> LF_disk_at_center.log
OMP_NUM_THREADS=1 ../../fargo start leapfrog_kdk.yml 1> LF_disk_at_center.log 2> LF_disk_at_center.log
OMP_NUM_THREADS=1 ../../fargo start euler.yml 1> euler_disk_at_center.log 2> euler_disk_at_center.log
OMP_NUM_THREADS=1 ../../fargo start nbody.yml 1> star_at_center.log 2> star_at_center.log
#./calc_deviation.py
