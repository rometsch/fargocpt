#!/usr/bin/env python3

import sys
sys.path.append("../../bin")
from fargocpt import run_fargo
import numpy as np
import re

Nproc_min = 1
Nproc_max = 1
Nthread_min = 1
Nthread_max = 2
for nproc in range(Nproc_min, Nproc_max+1):
    for nthread in range(Nthread_min, Nthread_max+1):

        filename = f"log_np{nproc}_nt{nthread}.txt" 
        with open(filename, "w") as logfile:
            run_fargo(nproc,nthread,["start", "testconfig.yml"], stdout=logfile, stderr=logfile)

        time_per_step_ms = np.genfromtxt("output/out/monitor/timestepLogging.dat", usecols=6)[1:]
        print(nproc, nthread, np.median(time_per_step_ms))
