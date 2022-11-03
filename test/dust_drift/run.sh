#!/usr/bin/env bash
make -C ../../src -j
if [[ $? > 0 ]]; then
    echo "Error in compilation. Refusing to continue!"
    exit $?
fi
# ../../fargo start dust_drift.yml
rm -rf output/dust_drift
mkdir -p output
sbatch -n 6 -o "output/%j.log" slurm_start.sh