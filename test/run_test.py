#!/usr/bin/env python3
import subprocess
import sys
import os
from argparse import ArgumentParser
import yaml

def compile_fargo(fargo_path, Nthreads, silent=False):

    if silent:
        out = open("make.log", "w")
        err = out
    else:
        out = sys.stdout
        err = sys.stderr

    cmd = ['make', '-j', str(Nthreads), '-C ' 'src/']
    subprocess.run(cmd, cwd=fargo_path, stdout=out, stderr=err)

    if silent: 
        out.close()

def run(fargo_path, setupfile, Nthreads, Nprocs, silent=False):
    import sys
    fargo_path = os.path.abspath(fargo_path)
    bindir = os.path.join(fargo_path, 'bin')
    sys.path.append(bindir)
    from fargocpt import run_fargo

    runname = os.path.basename(setupfile)[:-4]


    if silent:
        out = open(f"{runname}.out", "w")
        err = open(f"{runname}.err", "w")
    else:
        out = sys.stdout
        err = sys.stderr

    run_fargo(Nprocs, Nthreads, fargo_args=["start", setupfile], stdout=out, stderr=err)

    if silent:
        out.close()
        err.close()

def main():

    parser = ArgumentParser()
    parser.add_argument('-nt', type=int, default=2, help='Number of threads per process.')
    parser.add_argument('-np', type=int, default=1, help='Number of MPI processes.')
    parser.add_argument('-i', '--interactive', action='store_true', help='Interactive mode. Show the plot in a window.')
    parser.add_argument('-t', '--testonly', action='store_true', help='Only run the test.')
    parser.add_argument('-s', '--silent', action='store_true', help="Don't print anything.")

    opts = parser.parse_args()

    Nthreads = opts.nt
    Nprocs = opts.np

    with open("testconfig.yml", "r") as f:
        testconfig = yaml.safe_load(f)

    testname = testconfig["testname"]
    try:
        setupfiles = [testconfig["setupfile"]]
    except KeyError:
        setupfiles = testconfig["setupfiles"]


    if not opts.testonly:
        compile_fargo('../../', Nprocs*Nthreads, silent=opts.silent)
        for setupfile in setupfiles:
            run('../../', setupfile, Nthreads, Nprocs, silent=opts.silent)

    sys.path = [os.getcwd()] + sys.path
    from check_results import test
    test(f'../../output/tests/{testname}/out/')

if __name__=='__main__':
    main()