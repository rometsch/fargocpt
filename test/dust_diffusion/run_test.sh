#!/usr/bin/env bash
make -C src -j
mpirun -n 3 fargo start dust_test.par