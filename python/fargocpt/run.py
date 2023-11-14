#!/usr/bin/env python3

from argparse import ArgumentParser, RawTextHelpFormatter
from subprocess import Popen, PIPE, run
import signal
import sys
from sys import platform
from time import sleep
import os
import re
import tempfile
import yaml
import shutil

file_dir = os.path.abspath(os.path.dirname(__file__))
config_dir = os.path.expanduser("~/.config/fargocpt")
config_file = os.path.join(config_dir, "config.yml")


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
|-N <N> |                Perform N hydro steps


Following are the additional options available through this wrapper, refered to as [wrapper options] in the usage string:
"""

usage = "%(prog)s [wrapper options] [fargo options] start|restart <N>|auto configfile"

def main(args=sys.argv[1:]):
    opts, fargo_args = parse_opts(args)

    if opts.print_numa:
        print(get_numa_nodes_linux())
        print("Exiting now due to --print-numa option.")
        exit(0)

    try:
        executable_path = find_executable(opts.exe)
    except RuntimeError as e:
        print(e)
        exit(1)

    rv = run(fargo_args, np=opts.np, nt=opts.nt, mpi_verbose=opts.mpi_verbose,
              fallback_mpi=opts.fallback_mpi, fallback_openmp=opts.fallback_openmp,
              detach=opts.detach, exe=executable_path)
    sys.exit(rv)


def determine_cpu_allocation(np=None, nt=None):
    if np is not None and nt is not None:
        N_procs = np
        N_OMP_threads = nt
    elif np is None and nt is not None:
        N_procs = 1
        N_OMP_threads = nt
    elif np is not None and nt is None:
        N_procs = np
        N_OMP_threads = 1
    else:
        N_procs, N_OMP_threads = get_auto_num_procs_and_threads()
    return N_procs, N_OMP_threads

def detach_processGroup():
    """
    Detach a newly spawned process from the spawning process by
    setting the process group id to the process id.
    This will detach it from automatic signal propagation
    and allows fargo to exit gracefully.
    """
    os.setpgrp()

def find_executable(exe=None):
    """ Attempt to find the fargocpt executable.
    Check the following locations:
    1. Command line argument option
    2. ~/.config/fargocpt/config.yml
    3. In the path
    """
    executable_path = None
    if exe is not None:
        # 1. Command line argument option
        executable_path = exe
    
    if executable_path is None:
        # 2. ~/.config/fargocpt/config.yml
        try:
            with open(config_file, "r") as infile:
                config = yaml.safe_load(infile)
                executable_path = config["exe_path"]
        except FileNotFoundError:
            executable_path = None
        
    if executable_path is None:
        # 3. In the path
        guess = shutil.which("fargocpt_exe")
        if guess is not None:
            executable_path = guess

    if executable_path is None:
        msg = "Could not locate the fargocpt_exe executable.\nThere are three options to specify the path to it wich are checked in the following order.\n  1. specify it via the --exe cli option.\n  2. Create a config file at ~/.config/fargocpt/config.yml with the key 'exe_path' pointing to the fargocpt executable.\n  3. Put the bin/fargocpt_exe executable in your path and name it 'fargocpt_exe'."
        raise RuntimeError(msg)

    if not os.path.exists(executable_path):
        msg = f"Executable not found! I looked for {executable_path}\nPlease compile the code first: make -C <source dir in repo> -j 4"
        raise RuntimeError(msg)
    
    return executable_path


def parse_opts(args=sys.argv[1:]):
    parser = ArgumentParser(
        usage=usage, description=fargo_help_string, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-np", type=int, help="Number of processes to start.")
    parser.add_argument(
        "-nt", type=int, help="Number of threads to start per process.")
    parser.add_argument("--print-numa", action="store_true",
                        help="Print information about the numa nodes and exit.")
    parser.add_argument("--mpi-verbose", action="store_true",
                        help="Ask openmp and mpirun to output additional information about cpu allocation.")
    parser.add_argument("--fallback-mpi", action="store_true",
                        help="Fall back to a simpler call of mpi.")
    parser.add_argument("--fallback-openmp", action="store_true",
                        help="Fall back to calling openmp directly.")
    parser.add_argument("-d", "--detach", action="store_true", help="Detach the launcher from the simulation and run in the background.")
    parser.add_argument("-e", "--exe", type=str, default=None, help="Path to the fargocpt executable.")
    opts, remainder = parser.parse_known_args(args)
    return opts, remainder


def run(fargo_args, np=None, nt=None, mpi_verbose=False, stdout=None, stderr=None, fallback_mpi=False, fallback_openmp=False, envfile=None, detach=False, exe=None):
    """Run a fargo simulation.

    Args:
        np (int): Number of MPI processes to start.
        nt (int): Number of OpenMP threads to start per process.
        fargo_args (list of str): Arguments to be passed to the fargo executable.
        mpi_verbose (bool, optional): Verbose output of mpirun about process allocation. Defaults to False.
        stdout (file, optional): File to redirect stdout to. Defaults to None.
        stderr (file, optional): File to redirect stderr to. Defaults to None.
        fallback_mpi (bool, optional): If True, fall back to a simpler call of mpi. Defaults to False.
        fallback_openmp (bool, optional): If True, fall back to calling openmp directly. Defaults to False.
        envfile (str, optional): Before running fargo, source this file first. Use for loading modules on clusters.
        detach (bool, optional): Detach the launcher from the simulation and run in the background. Defaults to False.
        exe (str, optional): Path to the fargocpt executable. Defaults to None.
    """

    executable_path = find_executable(exe)
    N_procs, N_OMP_threads = determine_cpu_allocation(np=np, nt=nt)

    print(
        f"Running fargo with {N_procs} processes with {N_OMP_threads} OMP threads each using executable '{executable_path}'", flush=True)


    pidfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
    pidfilepath = pidfile.name
    pidfile.close()


    if not fallback_mpi and not fallback_openmp:
        cmd = ["mpirun"]
        cmd += ["-np", f"{N_procs}"]
        # write out pid of mpirun to correctly terminate multiprocess sims
        # if N_procs > 1:
        cmd += ["--report-pid", pidfilepath]
        if mpi_verbose:
            cmd += ["--display-map"]
            cmd += ["--display-allocation"]
            cmd += ["--report-bindings"]
        if mpi_verbose:
            cmd += ["-x", "OMP_DISPLAY_ENV=VERBOSE"]
        if N_OMP_threads > 1:
            if platform != "darwin":
                if os.path.exists("/.dockerenv"):
                    # bind to l3cache because bind to numa does not work inside docker
                    cmd += ["-x", "OMPI_MCA_rmaps_base_mapping_policy=l3cache"]
                    cmd += ["-x", "OMPI_MCA_hwloc_base_binding_policy=l3cache"]
                    cmd += ["--map-by", "ppr:1:socket"]
                    cmd += ["--bind-to", "socket"]
                else:
                    cmd += ["--map-by", "ppr:1:numa"]
                    cmd += ["--bind-to", "numa"]
            cmd += ["-x", "OMP_WAIT_POLICY=active"]
            cmd += ["-x", "OMP_PROC_BIND=close"]
            cmd += ["-x", "OMP_PLACES=cores"]
        cmd += ["-x", f"OMP_NUM_THREADS={N_OMP_threads}"]
        env = os.environ.copy()
        env_update = {}
    elif fallback_mpi:
        N_procs = N_procs * N_OMP_threads
        N_OMP_threads = 1
        cmd = ["mpirun"]
        cmd += ["--report-pid", pidfilepath]
        cmd += ["-np", f"{N_procs}"]
        env = os.environ.copy()
        env_update = {"OMP_NUM_THREADS": f"{N_OMP_threads}"}
        env.update(env_update)
    elif fallback_openmp:
        N_OMP_threads = N_procs*N_OMP_threads
        N_procs = 1
        env = os.environ.copy()
        env_update = {"OMP_NUM_THREADS": f"{N_OMP_threads}"}
        env.update(env_update)
        cmd = []
    cmd += [executable_path] + fargo_args
    if fallback_openmp:
        cmd += ["--pidfile", pidfilepath]
    
    if detach:
        if stdout is None:
            stdout = sys.stdout

        if stderr is None:
            stderr = sys.stderr

    cmd_desc = f"Running command: {' '.join(cmd)}"
    if len(env_update) > 0:
        cmd_desc += f" with env updated with {env_update}"
    print_wrapper(stdout, cmd_desc, flush=True)
    
    if envfile is not None:
        cmd = f"source {envfile}; " + " ".join(cmd)
    else:
        cmd = " ".join(cmd)

    p = Popen(cmd,
              stdout=stdout,
              stderr=stderr,
              env=env,
              preexec_fn=detach_processGroup,
              shell=True, executable="/bin/bash")

    # Read the pid from file
    for i in range(10):
        with open(pidfilepath, "r") as pidfile:
            pid = pidfile.read()
        if pid == "":
            sleep(0.1)
        else:
            break
    
    if not detach:
        print_wrapper(stdout, "fargo process pid", pid, flush=True)

        def handle_termination_request(signum, frame):
            try:
                import psutil
                pfargo = psutil.Process(int(pid))
                pfargo.send_signal(signal.SIGTERM)
            except ImportError:
                run(["kill", "-SIGTERM", f"{int(pid)}"])

        signal.signal(signal.SIGINT, handle_termination_request)
        signal.signal(signal.SIGTERM, handle_termination_request)

        if hasattr(p.stdout, "readline"):
            print("reading lines")
            for line in iter(p.stdout.readline, b''):
                line = line.decode("utf8")
                if stdout is not None and stdout != PIPE:
                    print_wrapper(stdout, line.rstrip())
                
        else:
            p.wait()

    else:
        print_wrapper(stdout, "fargo process pid", pid)
        print_wrapper(stdout, "detaching... check output dir for logs")

    return p.returncode

def print_wrapper(stdout, *args, **kwargs):
    if stdout != PIPE:
        print(*args, file=stdout, **kwargs)
    else:
        print(*args, **kwargs)

def get_num_cores():
    rv = None
    if "PBS_NP" in os.environ:
        try:
            rv = int(os.environ["PBS_NP"])
            print(f"Found PBS environment with {rv} cores")
        except KeyError:
            pass
    elif "SLURM_JOB_CPUS_PER_NODE" in os.environ and "SLURM_NNODES" in os.environ:
        try:

            # slurm_cpus_per_node = '72(x2),36' # example SLURM_JOB_CPUS_PER_NODE string for 2 nodes 72 cores each and 1 node with 36 cores
            slurm_cpus_per_node = os.environ["SLURM_JOB_CPUS_PER_NODE"]
            # we ignore nodes with different numbers of cores
            slurm_cpus_per_node = slurm_cpus_per_node.split(',')[0]
            # it is up to the user to assure that this does not cause problem (by only using nodes with equal number of cores)
            # filter out number of nodes used, this info is contained elsewhere
            slurm_cpus_per_node = slurm_cpus_per_node.split('(')[0]
            slurm_cpus_per_node = int(slurm_cpus_per_node)

            rv = int(os.environ["SLURM_NNODES"]) * slurm_cpus_per_node
            print(f"Found SLURM environment with {rv} cores")
        except KeyError:
            pass
    if rv is None:
        try:
            import psutil
            rv = psutil.cpu_count(logical=False)
        except ImportError:
            # Try to infer the number of cpus from the numa topology
            try:
                numa_nodes = get_numa_nodes_linux()
                rv = sum([len(v) for v in numa_nodes.values()])
                if rv == 0:
                    raise RuntimeError("Something went wrong while infering the number of available cpus. It shoud be > 0. Set it manually or install the 'psutil' python package!")
            except (FileNotFoundError):
                raise NotImplementedError("Infering number of available cpus from numa topology only works on linux systems. Set the number manually or install the 'psutil' python package.")
    return rv

def get_numa_nodes_linux():
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
    if not os.path.exists(path):
        raise NotImplementedError("This function only works on linux systems.")

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


def get_auto_num_procs_and_threads():
    """
    Determine a configuration that respects 
    the numa topology of the system.

    Returns:
        tuple: Number of processes and number of threads.
    """
    ncpu_avail = get_num_cores()

    try:
        numa_nodes = get_numa_nodes_linux()
        first_numa_node = [k for k in numa_nodes.keys()][0]
        N_cores_per_numa = len(numa_nodes[first_numa_node])
        N_procs = max(1, ncpu_avail//N_cores_per_numa)
        return N_procs, N_cores_per_numa
    except (NotImplementedError, FileNotFoundError, KeyError):
        return 1, ncpu_avail


if __name__ == "__main__":
    main()
