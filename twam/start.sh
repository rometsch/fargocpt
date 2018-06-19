#!/usr/bin/env bash

# BINAC: Change dir and load cuda
if [[ ${PBS_O_WORKDIR} && $(hostname -s) = gpu* ]]
then
	module load devel/cuda/8.0
	unset CUDA_VISIBLE_DEVICES
	cd $PBS_O_WORKDIR
else
	cd $PBS_O_WORKDIR
fi

src/fargo3d in/in.par
