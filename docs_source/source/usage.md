# Usage

How to use the code?

## Parallel execution

The code is parallelized using a hybrid MPI/OpenMP approach.
This means that the user has to specify the configuration of how the cores are distributed among the MPI processes and OpenMP threads.
The total number of processes x number of threads must be equal to the number of cores available on the machine.

In general, each compute node (physical computer/cpu socket/numa node) should have its own MPI process with a number of OpenMP threads equal to the number of cores on that node.

The `run_fargo` wrapper script takes the number of MPI proceeses via the `-np` argument and the number of OpenMP threads via the `-nt` argument.
Internally, `mpirun` is called with options suitable for an x86 processor architecture with numa nodes.

If no arguments are given, the wrapper script will try to determine the number of cores available on the machine and launch one MPI process per numa node with as many OpenMP threads as there are cores on that node.
The number of available cores is determined by in the following order:
1. Check whether there are environment variables set by the `SLURM` or `PBS` job schedulers.
2. Check the number of physical cores (no hyperthreading) using `psutil.cpu_count(logical=False)`.

If this does not work for you, there are also two fallback options which are less restrictive:
1. The `--fallback-mpi` will ignore OpenMP and just launch as many MPI processes as there are cores available (Nprocs * Nthreads, or maximum number of cores).
2. The `--fallback-openmp` will ignore MPI and just launch 1 process with as many OpenMP threads as there are cores available (Nprocs * Nthreads, or maximum number of cores).