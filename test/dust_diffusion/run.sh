#!/usr/bin/env bash
make -C ../../src -j
if [[ $? > 0 ]]; then
    echo "Error in compilation. Refusing to continue!"
    exit $?
fi

# oldpid="$(ps ax | grep fargo | grep mpirun | cut -d' ' -f1)"
# if [[ -n "$oldpid" ]]; then
#     kill $oldpid
# fi


LOGFILE="output/simulation.log"

rm -rf output/dust_diffusion
mkdir -p output

# if [[ -n "$(command -v squeue)" ]]
# then
#     sbatch -n 6 -o "$LOGFILE" slurm_start.sh
#     touch $LOGFILE
#     tail -f $LOGFILE
# else
#     ../../fargo start setup/simple_model.yml
# fi
export N_MPI=1
export OMP_NUM_THREADS=1
numactl mpirun -n $N_MPI --map-by NUMA:PE=$OMP_NUM_THREADS ../../fargo -v start setup/Charnoz2011_fast.yml
