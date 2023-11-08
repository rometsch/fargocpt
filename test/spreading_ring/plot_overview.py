#!/usr/bin/env python3

from argparse import ArgumentParser

import numpy as np
from scipy.special import iv
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument("outdir", type=str, help="Output directory to load data from.")
parser.add_argument("plotfilename", type=str, help="Filename of the plot.")
opts = parser.parse_args()

def calc_theo(r, t):
    R0 = 1
    M = 1.0
    nu = 4.77e-5 / R0**2

    x = r / R0
    tau0 = 0.016
    tau = 12 * nu * t / R0**2 + tau0

    I = iv(0.25, 2.0*x/tau)
    Sigma = M / (np.pi * R0**2) / tau / x**(1/4) * I * np.exp(-(1+x**2)/tau)
    I0 = iv(0.25, 2.0*x/tau0)
    Sigma0 = M / (np.pi * R0**2) / tau0 / x**(1/4) * I0 * np.exp(-(1+x**2)/tau0)
    
    return Sigma, Sigma0


dt = 1
out = opts.outdir

fig, axs = plt.subplots(2,1,figsize=(6,4.5))
axs = np.ravel(axs)

Radii = np.loadtxt(out + "/used_rad.dat", skiprows=0)
Radii = Radii
Rinf = Radii[:-1]
Rsup = Radii[1:]
Rmed = 2.0/3.0*(Rsup*Rsup*Rsup-Rinf*Rinf*Rinf)
Rmed = Rmed / (Rsup*Rsup-Rinf*Rinf)

Quantities = np.loadtxt(out + '/monitor/Quantities.dat', skiprows=26)
t = Quantities[int(dt),2]
Sigma, Sigma0 = calc_theo(Rmed, t)

axs[0].plot(Rmed, Sigma, ls='-', color='C0', lw=2.5, label='Theory')
axs[0].plot(Rmed, Sigma0, ls='--', color='black', lw=2.5, label='Initial')

file_name = out + f'/snapshots/{dt}/Sigma.dat'
data = np.fromfile(file_name)
N = len(data)
nr = len(Rmed)
nphi = int(N/nr)

data = data.reshape((nr, nphi))

data = np.mean(data, 1)

axs[0].plot(Rmed, data, ls='--', color='C1', lw=1.5, label='Simulation')

diff = np.abs(Sigma / data - 1)

axs[1].plot(Rmed, diff, ls='-', color='C0', lw=2.5)

axs[0].legend(loc='upper right')

ax = axs[0]
ax.set_ylabel('$\Sigma$ [code units]')
ax.set_xlabel('x')

ax = axs[1]
ax.set_ylabel('$\Sigma$ Relative difference')
ax.set_xlabel('x')
ax.set_yscale('log')

plt.savefig(opts.plotfilename, dpi=150, bbox_inches='tight')