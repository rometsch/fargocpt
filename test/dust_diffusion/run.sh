#!/usr/bin/env bash
make -C ../../src -j
if [[ $? > 0 ]]; then
    echo "Error in compilation. Refusing to continue!"
    exit $?
fi

rm -rf output/dust_diffusion
mkdir -p output

if [[ -n "$(command -v squeue)" ]]
then
    sbatch -n 6 -o "output/simulation.log" slurm_start.sh
    tail -f output/simulation.log
else
    ../../fargo start setup/Charnoz2011_fast.yml
fi