#!/usr/bin/env bash

echo "Prelaunch config for fargocpt"

NCPUS_PER_NODE=$(lscpu | grep "^\s*CPU(s):" | tr " " "\n" | tail -n 1)
NTHREADS_PER_CORE=$(lscpu | grep "^\s*Thread(s) per core:" | tr " " "\n" | tail -n 1)
NNUMA_NODES=$(lscpu | grep "NUMA node(s):" | tr " " "\n" | tail -n 1)
N_THREADS_PER_NUMA=$(python3 -c "print($NCPUS_PER_NODE//$NNUMA_NODES//$NTHREADS_PER_CORE)")

echo "Detected $NCPUS_PER_NODE cpus per node"
echo "Detected $NTHREADS_PER_CORE threads per core"
echo "Detected $NNUMA_NODES numa nodes"
echo "Using $N_THREADS_PER_NUMA OMP threads per process"


if [[ -n ${PBS_NP+x} ]];
then
	N_CORES_AVAIL=$PBS_NP
	echo "Detected PBS environment with $PBS_NP cores available."
	echo "Running on hosts $(cat $PBS_NODEFILE | uniq | sort -V | tr '\n' ', ')"
elif [[ -n ${SLURM_NPROCS+x} ]];
then
	N_CORES_AVAIL=$SLURM_NPROCS
	echo "Detected SLURM environment with $SLURM_NPROCS cores available."
	echo "Running on hosts $SLURM_NODELIST"
else
	N_CORES_AVAIL=$(python3 -c "print($NCPUS_PER_NODE//$NTHREADS_PER_CORE)")
	echo "Did not detect PBS or SLURM queuing system. I'll be greedy and take all $N_CORES_AVAIL systems cores..."
fi

if (( $N_CORES_AVAIL < $N_THREADS_PER_NUMA ));
then
	NPROCS=1
	N_OMP_THREADS=$N_CORES_AVAIL
else
	NPROCS=$(python3 -c "print($N_CORES_AVAIL//$N_THREADS_PER_NUMA)")
	N_OMP_THREADS=$N_THREADS_PER_NUMA
fi


echo "Starting $NPROCS processes with $N_THREADS_PER_NUMA OMP Threads, current host $(hostname)"

STARTCMD="./fargo $@"

mpirun \
	--np $NPROCS \
	--display-map \
	--display-allocation \
	--report-bindings \
	--map-by ppr:1:numa \
	--bind-to numa \
	-x OMP_DISPLAY_ENV=VERBOSE \
	-x OMP_WAIT_POLICY=active \
	-x OMP_PROC_BIND=close \
	-x OMP_PLACES=cores \
	-x OMP_NUM_THREADS=$N_OMP_THREADS \
	$STARTCMD
