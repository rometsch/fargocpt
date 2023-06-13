#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
../../run_fargo -nt 2 -np 1 start setup.yml 1> out.log 2>err.log
./calc_deviation.py
