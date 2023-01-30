#!/usr/bin/env python3

from argparse import ArgumentParser, RawTextHelpFormatter
from subprocess import run
import os
import re

file_dir = os.path.abspath(os.path.dirname(__file__))
executable_path = os.path.join(file_dir, "fargocpt")


fargo_help_string = """This script is a wrapper to start the fargocpt code.
Its main use is to specify the number of processes started and threads used per process, and to automatically infer some information about the available cpus and the numa nodes on the system.
Beware: this wrapper is only tested with open-mpi. Other flavors might fail.

Help string from the fargo executable:
| 
| start                  Start a new simulation from scratch
| restart <N>            Restart from an old simulation output, latest if no N specified
| auto                   Same as restart if output files are present, otherwise same as start
| -d | --debug           Print some debugging information on 'stdout' at each timestep
| -v | --verbose         Verbose mode. Tells everything about parameters file
| -q | --quiet           Only print errors and warnings
| -b |                   Adjust azimuthal velocity to impose strict centrifugal balance at t=0
| -c |                   Sloppy CFL condition (checked at each DT, not at each timestep)
| -n |                   Disable simulation. The program just reads parameters file
| -m |                   estimate memory usage and print out


Following are the additional options available through this wrapper, refered to as [wrapper options] in the usage string:
"""

usage="%(prog)s [wrapper options] [fargo options] start|restart <N>|auto configfile"

def main():
    opts, fargo_args = parse_opts()

    if opts.print_numa:
        print(get_numa_nodes())
        print("Exiting now due to --print-numa option.")
        exit(0)
        
    if not os.path.exists(executable_path):
        print(f"Executable not found! I looked for {executable_path}.")
        src_dir = os.path.relpath(os.path.join(
            file_dir, "..", "src"), os.getcwd())
        print(f"Please compile the code first: make -C {src_dir} -j 4")
        exit(1)

    ncpu = get_num_cores()
    numa_nodes = get_numa_nodes()
    first_numa_node = [k for k in numa_nodes.keys()][0]
    N_cores_per_numa = len(numa_nodes[first_numa_node])

    if opts.np is None and opts.nt is None:        
        N_procs = max(1, ncpu//N_cores_per_numa)
        N_OMP_threads = ncpu//N_procs

    elif opts.np is not None and opts.nt is None:
        N_procs = opts.np
        N_OMP_threads = N_cores_per_numa
    elif opts.np is None and opts.nt is not None:
        N_procs = 1
        N_OMP_threads = opts.nt
    else:
        N_procs = opts.np
        N_OMP_threads = opts.nt

    print(
        f"Running fargo with {N_procs} procs with {N_OMP_threads} OMP threads each.")

    run_fargo(N_procs, N_OMP_threads, fargo_args, mpi_verbose=opts.mpi_verbose)


def parse_opts():
    parser = ArgumentParser(usage=usage, description=fargo_help_string, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-np", type=int, help="Number of processes to start.")
    parser.add_argument(
        "-nt", type=int, help="Number of threads to start per process.")
    parser.add_argument("--print-numa", action="store_true",
                        help="Print information about the numa nodes and exit.")
    parser.add_argument("--mpi-verbose", action="store_true",
                        help="Ask openmp and mpirun to output additional information about cpu allocation.")
    opts, remainder = parser.parse_known_args()
    return opts, remainder


def run_fargo(N_procs, N_OMP_threads, fargo_args, mpi_verbose=False, stdout=None, stderr=None):
    """Run a fargo simulation.

    Args:
        N_procs (int): Number of MPI processes to start.
        N_OMP_threads (int): Number of OpenMP threads to start per process.
        fargo_args (list of str): Arguments to be passed to the fargo executable.
        mpi_verbose (bool, optional): Verbose output of mpirun about process allocation. Defaults to False.
    """
    cmd = ["mpirun"]
    cmd += ["--np", f"{N_procs}"]
    if mpi_verbose:
        cmd += ["--display-map"]
        cmd += ["--display-allocation"]
        cmd += ["--report-bindings"]
    if mpi_verbose:
        cmd += ["-x", "OMP_DISPLAY_ENV=VERBOSE"]
    if N_OMP_threads > 1:
        cmd += ["--map-by", "ppr:1:numa"]
        cmd += ["--bind-to", "numa"]
        cmd += ["-x", "OMP_WAIT_POLICY=active"]
        cmd += ["-x", "OMP_PROC_BIND=close"]
        cmd += ["-x", "OMP_PLACES=cores"]
    cmd += ["-x", f"OMP_NUM_THREADS={N_OMP_threads}"]
    cmd += [executable_path] + fargo_args
    run(cmd, stdout=stdout, stderr=stderr)


def get_num_cores():
    try:
        rv = int(os.environ["PBS_NP"])
        print(f"Found PBS environment with {rv} cores")
    except KeyError:
        try:
            rv = int(os.environ["SLURM_NNODES"]) * \
                int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
            print(f"Found SLURM environment with {rv} cores")
        except KeyError:
            rv = get_num_cores_from_numa()
            print(
                f"Found no PBS or SLURM environment.\nI'll be greedy and take all {rv} system cores!")
    return rv


def get_num_cores_from_numa():
    """Return the number of physical cores on the system.

    Hyperthreads are ignored, as they don't speed up the computation.

    Returns:
        int: Number of cores.
    """
    numa_nodes = get_numa_nodes()
    rv = sum([len(n) for n in numa_nodes.values()])
    return rv


def get_numa_nodes():
    """Return the numa topology of the system.

    Data is extracted from /sys/devices/system/cpu.
    The numa nodes on the system are the keys of the dictionary.
    Each value is a list of tuples which themselves represent physical cores.
    Each tuples contains the id of the threads belonging to the physical core.
    Without hyperthreading, each tuple contains one id, with hyperthreading each tuple contains two keys.

    Returns:
        dict: Topology of the numa nodes.
    """

    path = "/sys/devices/system/cpu"

    cpus = [d for d in os.listdir(path) if re.match("cpu\d+", d) is not None]

    nodes = {}

    for cpu in cpus:
        cpu_path = os.path.join(path, cpu)
        cpu_node = [s for s in os.listdir(
            cpu_path) if re.match("node\d+", s)][0]
        with open(os.path.join(cpu_path, "topology", "thread_siblings_list"), "r") as infile:
            threads = sorted(infile.read().strip().split(","))
        threads = tuple(threads)
        try:
            nodes[cpu_node].add(threads)
        except KeyError:
            nodes[cpu_node] = set()
            nodes[cpu_node].add(threads)

    for key in nodes:
        nodes[key] = sorted(nodes[key])

    return nodes


if __name__ == "__main__":
    main()
