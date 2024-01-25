#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import yaml

def main():

    with open("torque_test.yml", "r") as infile:
        params = yaml.safe_load(infile)

    alpha = params["SigmaSlope"] #1.5 # - surface density power law index
    flaringindex = params["FlaringIndex"]
    beta = 1-2*flaringindex #1 # - temperature power law index
    q = float(params["nbody"][1]["mass"])/float(params["nbody"][0]["mass"]) #2e-5 # reduced planet mass
    h = params["AspectRatio"]#0.05 # aspect ratio
    SigmaP = 3.76e-4 # Surface density at planet location
    OmegaP = 1#2*np.pi
    b = h*params["ThicknessSmoothing"] #0.4*h # smoothing lenght factor

    Gamma0 = (q/h)**2 * SigmaP * OmegaP**2
    # Paardekooper et al 2009

    expected_torque = -(2.5 + 1.7*beta - 0.1*alpha) * (0.4/(b/h))**0.71

    # print("alpha", alpha)
    # print("beta", beta)
    # print("q", q)
    # print("h", h)
    # print("SigmaP", SigmaP)
    # print("OmegaP", OmegaP)
    # print("b", b)
    # print("expected torque", expected_torque)


    fig, ax = plt.subplots(figsize=(6,4))

    dirs = []
    for name in os.listdir("../../output/tests/planet_torque"):
        dirs.append(f"../../output/tests/planet_torque/{name}")
        
    for outdir in dirs:
        # print(outdir)

        fname = f"{outdir}/monitor/nbody1.dat"
        name = outdir.rstrip("/").split("/")[-1]
        try:
            time, torque, a = np.genfromtxt(fname, usecols=(7,18,12), unpack=True)
        except OSError:
            print(f"{fname} not found")
            continue
        time = time/(2*np.pi)
        ax.plot(time, torque/Gamma0, label=name)

    ax.set_xlabel(r"$t$ [orbits]")
    ax.set_ylabel(r"$\Gamma / \Gamma_0$")


    ax.axhline(expected_torque, color="gray", label="theory")
    ax.legend()

    fig.savefig("plot.jpg", dpi=150, bbox_inches='tight')
    fig.savefig("plot.pdf", dpi=300, bbox_inches='tight')

    ratios = torque/Gamma0 / expected_torque
    average = np.average(ratios[-len(ratios)//10:])

    threshold = 0.2

    with open("test_log.txt", "w") as of:
        print(f"Deviation over last 10 percent of time is {average-1} and threshold is {threshold}.", file=of)

    test_name = os.path.basename(os.getcwd())
    if (np.abs(average - 1) < threshold):
        print(f"SUCCESS: {test_name}")
    else:
        print(f"FAIL: {test_name}")

    with open("test.log", "w") as f:
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{current_time}", file=f)
        print(f"Deviation over last 10 percent of time is {average-1} and threshold is {threshold}.", file=f)


if __name__=="__main__":
    main()
