#!/usr/bin/env python3
import subprocess
import sys
import os
from argparse import ArgumentParser
import yaml
from fargocpt import run as run_fargo

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

    runname = os.path.basename(setupfile)[:-4]

    if silent:
        out = open(f"{runname}.out", "w")
        err = open(f"{runname}.err", "w")
    else:
        out = sys.stdout
        err = sys.stderr

    returncode = run_fargo(fargo_args=["start", setupfile], np=Nprocs, nt=Nthreads, stdout=out, stderr=err, exe=f"{fargo_path}/bin/fargocpt_exe")

    if silent:
        out.close()
        err.close()

    return returncode

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


    errors_happened = False
    if not opts.testonly:
        compile_fargo('../../', Nprocs*Nthreads, silent=opts.silent)
        for setupfile in setupfiles:
            returncode = run('../../', setupfile, Nthreads, Nprocs, silent=opts.silent)
            if returncode != 0:
                print(f"ERROR: Run {setupfile} failed with return code {returncode}.")
                errors_happened = True

    if errors_happened:
        print("ERROR: Some runs failed.")
        sys.exit(1)

    sys.path = [os.getcwd()] + sys.path
    from check_results import test
    test(f'../../output/tests/{testname}/out/')

if __name__=='__main__':
    main()