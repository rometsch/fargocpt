#!/usr/bin/env bash

# sbatch --sockets-per-node=2 --cores-per-socket=6 slurm_start.sh
mkdir -p output
rv=$(sbatch -n 6 -o "output/%j.log" slurm_start.sh)
number=$(echo $rv | awk '{print $4}')
sleep 2
tail -f slurm-$number.out
