#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as units

def main():

    alpha = 1.5 # - surface density power law index
    beta = 1 # - temperature power law index
    q = 2e-5 # reduced planet mass
    h = 0.05 # aspect ratio
    SigmaP = 3.76e-4 # Surface density at planet location
    OmegaP = 1#2*np.pi
    b = 0.4*h # smoothing lenght factor

    Gamma0 = (q/h)**2 * SigmaP * OmegaP**2
    # Paardekooper et al 2009

    expected_torque = -(2.5 + 1.7*beta - 0.1*alpha) * (0.4/(b/h))**0.71

    fig, ax = plt.subplots(figsize=(10,8))

    fname = "output/monitor/bigplanet2.dat"
    time, torque = np.genfromtxt(fname, usecols=(6,17), unpack=True)
    time = time/(2*np.pi)
    ax.plot(time, torque/Gamma0, label="corr")

    # fname = "output-nosgcorr/" + "bigplanet2.dat"
    # x, y, time, torque = np.genfromtxt(fname, usecols=(1,2,6,17)).T
    # time = time/(2*np.pi)
    # ax.plot(time, torque/Gamma0, label="no-corr")


    ax.axhline(expected_torque, color="C2", label="theory")
    ax.set_xlabel(r"$t$ [orbits]")
    ax.set_ylabel(r"$\Gamma / \Gamma_0$")

    ax.legend()

    plt.show()

if __name__=="__main__":
    main()