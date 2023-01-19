#!/usr/bin/env python3

import sys
sys.path.append("../../bin")
from fargocpt import run_fargo
import numpy as np
import re

Nproc_min = 1
Nproc_max = 2
Nthread_min = 1
Nthread_max = 6
for nproc in range(Nproc_min, Nproc_max+1):
    for nthread in range(Nthread_min, Nthread_max+1):

        filename = f"log_np{nproc}_nt{nthread}.txt" 
        with open(filename, "w") as logfile:
            run_fargo(nproc,nthread,["start", "testconfig.yml"], stdout=logfile, stderr=logfile)


        with open(filename, "r") as logfile:
            ptrn = re.compile(r"Time per Step:\s([0-9\.]+)\smilliseconds")
            for line in logfile:
                m = re.search(ptrn, line)
                if m is not None:
                    ms_per_step = m.groups()[0]

        # num_steps, time = np.genfromtxt("output/out/monitor/timestepLogging.dat", usecols=(4,3))[1:].T
        # time_per_hydrostep = (time[-1] - time[0])/np.sum(num_steps)
        print(nproc, nthread, ms_per_step)
