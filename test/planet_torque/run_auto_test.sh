#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR
../../run_fargo -nt 2 -np 1 start -N 1000 torque_test.yml 1> out.log 2>err.log
../../run_fargo -nt 2 -np 1 restart torque_test.yml 1> out.log 2>err.log
./plot_torque.py
