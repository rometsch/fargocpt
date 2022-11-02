#!/usr/bin/env bash
rm ../../fargo
make -C ../../src -j
# ../../fargo start dust_drift.yml
sbatch -n 6 -o "output/%j.log" slurm_start.sh