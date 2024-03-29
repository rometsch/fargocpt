#!/usr/bin/env python3

import os
import numpy as np
from types import SimpleNamespace
from argparse import ArgumentParser
from matplotlib import colormaps

import yaml

def main():

    test_name = "FLD1D"
    success = calc_deviation("../../output/tests/FLD1D/out")

    if success:
        print(f"SUCCESS: {test_name}")
    else:
        print(f"FAIL: {test_name} at {os.getcwd()}")




def calc_deviation(outdir):

    data = load_data(outdir)

    Nfirst = data.Ns[0]
    Nlast = data.Ns[-1]


    Rmin, Rmax = 0.2, 10
    R = np.geomspace(Rmin, Rmax, 1000)
    theo = theoretical_results(data.rc)

    deltaT = data.Tprofiles[Nlast] / theo.T - 1
    dev = np.max(np.abs(deltaT[data.rc < 9.5]))

    success = dev < 0.1

    with open("test.log", "w") as f:
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{current_time}", file=f)
        print(f"max deviation = {dev}, theshold = 0.1", file=f)

    return success


def theoretical_results(R):
    """Calculate the theoretical equilibrium temperature profile.

    Credit: Alex Ziampras

    Args:
        R: Radius array in units of au.

    Returns:
        Namespace with results
    """
    mu = 2.353
    K = 106701.29 # code unit for temperature
    T0 = mu * 0.05**2 * K

    Rmin = R[0]
    Rmax = R[-1]
    
    f1, f2 = -3.5, 5 # old module
    f1, f2 = -2, 9/2 # new module

    R1 = Rmin ** f1
    R2 = Rmax ** f1
    T1 = (T0 / Rmin)**f2
    T2 = (T0 / Rmax)**f2
    c1 = (T2-T1) / (R2-R1)
    c2 = (R2*T1 - R1*T2) / (R2-R1)
    T = (c1 * R ** f1 + c2) ** (1/f2)

    from types import SimpleNamespace

    theo = SimpleNamespace()
    theo.R = R
    theo.T = T
    theo.T0 = T0
    theo.T0p = T0/R
    theo.Sigma = R**-0.5
    theo.H = T**0.5*R**1.5
    return theo


def load_data(outdir):

    data = SimpleNamespace()

    data.Nr, data.Naz = np.genfromtxt(
        f"{outdir}/dimensions.dat", usecols=(4, 5), unpack=True, dtype=int)
    data.ri = np.genfromtxt(f"{outdir}/used_rad.dat")  # [1:-1]
    data.phii = np.linspace(0, 2*np.pi, data.Naz+1)
    data.rc = 0.5*(data.ri[1:] + data.ri[:-1])

    data.Ns = np.genfromtxt(f"{outdir}/snapshots/list.txt", dtype=int)

    tempunit = 1
    with open(f"{outdir}/units.yml", "r") as infile:
        units = yaml.safe_load(infile)
        tempunit = units["temperature"]["cgs value"]

    data.Ts = {n: tempunit*np.fromfile(f"{outdir}/snapshots/{n}/Temperature.dat",
                                       dtype=np.float64).reshape(data.Nr, data.Naz) for n in data.Ns}
    data.Tprofiles = {n: np.average(T, axis=1) for n, T in data.Ts.items()}

    data.Sigmas = {n: np.fromfile(f"{outdir}/snapshots/{n}/Sigma.dat",
                                  dtype=np.float64).reshape(data.Nr, data.Naz) for n in data.Ns}
    data.Sigmaprofiles = {n: np.average(S, axis=1)
                          for n, S in data.Sigmas.items()}

    return data

if __name__ == "__main__":
    main()
