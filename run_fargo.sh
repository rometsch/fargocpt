#!/usr/bin/env bash

echo "Prelaunch config for fargocpt"

NCPUS_PER_NODE=$(lscpu | grep "^\s*CPU(s):" | tr " " "\n" | tail -n 1)
NTHREADS_PER_CORE=$(lscpu | grep "^\s*Thread(s) per core:" | tr " " "\n" | tail -n 1)
NNUMA_NODES=$(lscpu | grep "NUMA node(s):" | tr " " "\n" | tail -n 1)
NOMP_THREADS=$(python3 -c "print($NCPUS_PER_NODE//$NNUMA_NODES//$NTHREADS_PER_CORE)")

echo "Detected $NCPUS_PER_NODE cpus per node"
echo "Detected $NTHREADS_PER_CORE threads per core"
echo "Detected $NNUMA_NODES numa nodes"
echo "Using $NOMP_THREADS OMP threads per process"


if [[ -n "${PBS_NP+x}" ]];
then
	N_CORS_AVAIL=$PBS_NP
	echo "Detected PBS environment with $PBS_NP cores available."
else
	N_CORES_AVAIL=$(python3 -c "print($NCPUS_PER_NODE//$NTHREADS_PER_CORE)")
	echo "Did not detect any queuing system. I'll be greedy and take all the systems cores..."
fi
NPROCS=$(python3 -c "print($N_CORES_AVAIL//$NOMP_THREADS)")

echo "Starting $NPROCS processes with $NOMP_THREADS OMP Threads, current host $(hostname)"

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
	-x OMP_NUM_THREADS=$NOMP_THREADS \
	$STARTCMD
