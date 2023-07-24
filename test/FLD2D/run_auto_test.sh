#!/usr/bin/env bash

FILEDIR="$(dirname $(realpath $0))"
cd $FILEDIR

./run_code.sh
./check_solution.py --diff