#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR

./run_code.sh &> simulation.log
./check_solution.py --test
./check_solution.py --outfile plot.jpg