#!/usr/bin/env python3
import os

testname = 'mpi_simple'

def test(outdir):
    identifier = os.path.join(outdir, "snapshots/1/misc.bin")
    if os.path.exists(identifier):
        print(f"SUCCESS: {testname}")
    else:
        print(f"FAILURE: {testname}")

