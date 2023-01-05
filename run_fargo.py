#!/usr/bin/env python3

from argparse import ArgumentParser
from subprocess import run
import os
import re
from multiprocessing import cpu_count

def main():
    opts = parse_opts()
    
    if opts.print_numa:
        print(get_numa_nodes())
    
    if opts.np is None and opts.nt is None:
        ncpu = get_num_cores()
        numa_nodes = get_numa_nodes()
        first_numa_node = [k for k in numa_nodes.keys()][0]
        N_cores_per_numa = len(numa_nodes[first_numa_node])
        N_procs = max(1, ncpu//N_cores_per_numa)
        N_OMP_threads = ncpu//N_procs
            
    elif opts.np is not None and opts.nt is None:
        N_procs = opts.np
        N_OMP_threads = 1
    elif opts.np is None and opts.nt is not None:
        N_procs = 1
        N_OMP_threads = opts.nt
    else:
        N_procs = opts.np
        N_OMP_threads = opts.nt

    print(f"Running fargo with {N_procs} procs with {N_OMP_threads} OMP threads each.")
    
    start(opts.configfile, opts.mode, N_procs, N_OMP_threads)

def parse_opts():
    parser = ArgumentParser()
    parser.add_argument("mode", choices=["start", "restart", "auto"], help="Start mode for fargo.")
    parser.add_argument("configfile", help="Path of the config file.")
    parser.add_argument("-np", type=int, help="Number of processes to start.")
    parser.add_argument("-nt", type=int, help="Number of threads to start per process.")
    parser.add_argument("--print-numa", action="store_true")
    opts = parser.parse_args()
    return opts

def start(configfile, mode, N_procs, N_OMP_threads):

    cmd = ["mpirun"]
    cmd += ["--np", f"{N_procs}"]
    cmd += ["--display-map"]
    cmd += ["--display-allocation"]
    cmd += ["--report-bindings"]
    cmd += ["--map-by", "ppr:1:numa"]
    cmd += ["--bind-to", "numa"]
    cmd += ["-x", "OMP_DISPLAY_ENV=VERBOSE"]
    cmd += ["-x", "OMP_WAIT_POLICY=active"]
    cmd += ["-x", "OMP_PROC_BIND=close"]
    cmd += ["-x", "OMP_PLACES=cores"]
    cmd += ["-x", f"OMP_NUM_THREADS={N_OMP_threads}"]
    cmd += ["./fargo", mode, configfile]
    run(cmd)

def get_num_cores():
    try:
        rv = os.environ["PBS_NP"]
        print(f"Found PBS environment with {rv} cores")
    except KeyError:
        try:
            rv = os.environ["SLURM_NPROCS"]
            print(f"Found SLURM environment with {rv} cores")
        except KeyError:
            rv = get_num_cores_from_numa()
            print(f"Found no PBS or SLURM environment.\nI'll be greedy and take all {rv} system cores!")
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
        cpu_node = [s for s in os.listdir(cpu_path) if re.match("node\d+",s)][0]
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

if __name__=="__main__":
    main()