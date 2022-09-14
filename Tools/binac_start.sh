#!/usr/bin/env bash

if [[ ${PBS_O_WORKDIR} ]]; then
    cd $PBS_O_WORKDIR
fi

N_MPI=4
N_OMPT=7

module load numlib/fftw/3.3.8-openmpi-3.1-gnu-9.2
module load mpi/openmpi/3.1-gnu-9.2
module load numlib/gsl/2.5-gnu-9.2

export OMPI_MCA_mpi_cuda_support=0
export OMP_DISPLAY_ENV=VERBOSE
export OMP_NUM_THREADS=$N_OMPT
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_WAIT_POLICY=active

numactl mpirun -n $N_MPI --map-by NUMA:PE=$OMP_NUM_THREADS ./fargo start testconfig.yml
