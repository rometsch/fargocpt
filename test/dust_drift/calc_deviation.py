#!/usr/bin/env python3

import os
import numpy as np

from drift_theo import vdrift_theo
from load_dust import construct_dust_trajectories

header = """Deviations of drift speeds in simulation to theoretical drift speeds (sim/theo -1)
for different stokes numbers averaged over the last tenth of the time series.
Syntax: Stokesnumber\tdeviation"""


def main():
    success, res = calc_deviation("../../output/tests/dust_drift/out")
    np.savetxt("deviations.txt", res, header=header)
    
    
    test_name = "Dust drift"

    if success:
        print(f"SUCCESS: {test_name}")
    else:
        print(f"FAIL: {test_name} at {os.getcwd()}")
        print(f"See {os.getcwd()}/deviations.txt for more detail.")

def calc_deviation(outdir):

    particles = construct_dust_trajectories(outdir)

    tolerance = 0.01

    successes = []
    
    res = np.zeros((len(particles), 2))
    
    for i, p in particles.items():
        t = p["time"]
        r = p["r"]

        stokes = p["stokes"]
        vtheo = vdrift_theo(stokes, r.to("au"))
        vtheo = vtheo.to("cm/s")
                    
        r = r.to("cm")
        t = t.to("s")
        rdot = (r[1:] - r[:-1])/(t[1:] - t[:-1])
        rdot = rdot.to("cm/s")
                    
        Navg = len(rdot)//10
        q = np.average(rdot[-Navg:]/vtheo[-Navg:])
        
        St = np.average(stokes[-Navg:])
        
        success = np.abs(q - 1) < tolerance
        
        res[i] = (St, q -1)
        
        successes.append(success)


    return all(successes), res

if __name__ == "__main__":
    main()
