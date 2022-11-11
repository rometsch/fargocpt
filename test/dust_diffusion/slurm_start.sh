#!/usr/bin/env bash

N_MPI=1
N_OMPT=6

# export OMPI_MCA_mpi_cuda_support=0
export OMP_DISPLAY_ENV=VERBOSE
export OMP_NUM_THREADS=$N_OMPT
# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
# export OMP_WAIT_POLICY=active

# numactl mpirun -n $N_MPI --map-by NUMA:PE=$OMP_NUM_THREADS ../../fargo -v start dust_drift.yml
numactl mpirun -n $N_MPI --map-by NUMA:PE=$OMP_NUM_THREADS ../../fargo -v start setup/simple_model.yml
