#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
../../run_fargo -nt 2 -np 1 start torque_test.yml 1> out.log 2>err.log
./plot_torque.py
