#!/usr/bin/env python3

import os

import numpy as np
from scipy.special import iv
from types import SimpleNamespace

def load_data(outdir):

    data = SimpleNamespace()

    data.Nr, data.Naz = np.genfromtxt(
        f"{outdir}/dimensions.dat", usecols=(4, 5), unpack=True, dtype=int)
    data.ri = np.genfromtxt(f"{outdir}/used_rad.dat")  # [1:-1]
    data.phii = np.linspace(0, 2*np.pi, data.Naz+1)
    
    Radii = data.ri
    Rinf = Radii[:-1]
    Rsup = Radii[1:]
    Rmed = 2.0/3.0*(Rsup*Rsup*Rsup-Rinf*Rinf*Rinf)
    Rmed = Rmed / (Rsup*Rsup-Rinf*Rinf)
    
    data.rc = Rmed

    data.Ns = np.genfromtxt(f"{outdir}/snapshots/list.txt", dtype=int)
    n = data.Ns[-1]

    Sigma = np.fromfile(f"{outdir}/snapshots/{n}/Sigma.dat", dtype=np.float64).reshape(data.Nr, data.Naz)
    data.Sigma = np.mean(Sigma, 1)
    
    dt = n
    Quantities = np.loadtxt(outdir + '/monitor/Quantities.dat', skiprows=26)
    data.t = Quantities[int(dt),2]

    return data

def calc_theo(r, t):
    R0 = 1
    M = 1.0
    nu = 4.77e-5 / R0**2

    x = r / R0
    tau0 = 0.016
    tau = 12 * nu * t / R0**2 + tau0

    I = iv(0.25, 2.0*x/tau)
    Sigma = M / (np.pi * R0**2) / tau / x**(1/4) * I * np.exp(-(1+x**2)/tau)
    
    return Sigma

def calc_deviation(outdir):

    data = load_data(outdir)
    theo = calc_theo(data.rc, data.t)
    
    diff = np.abs(data.Sigma / theo - 1)
    
    average = np.mean(diff)
 
    threshold = 0.0045
    success = np.abs(average) < threshold
    
    return success   
    
def main():

    success = calc_deviation("../../output/tests/spreading_ring/out")
    
    
    test_name = os.path.basename(os.getcwd())    
    if (success):
        print(f"SUCCESS: {test_name}")
    else:
        print(f"FAIL: {test_name}")

if __name__=='__main__':
    main()


