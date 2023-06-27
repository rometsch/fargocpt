#!/usr/bin/env python3

import numpy as np
import yaml
import astropy.constants as const
import astropy.units as units


# get size of grid specified in the setup file
with open("setup.yml", "r") as infile:
    params = yaml.safe_load(infile)
    Nrad = params["Nrad"]
    Naz = params["Nsec"]
    rmin = params["Rmin"]
    rmax = params["Rmax"]

ri = np.geomspace(rmin, rmax, Nrad+1)
phii = np.linspace(0, 2*np.pi, Naz+1)
Ri, Phii = np.meshgrid(ri, phii, indexing="ij")

dphi = phii[1] - phii[0]
A = 0.5*(Ri[1:,1:]**2 - Ri[:-1,1:]**2)*dphi

# now create an energy array with values in cgs with the energy located in one cell

nr = 3*Nrad//4
nphi = 0

T0 = 100*units.K
cv = const.ch
# energy0 = 

energy0 = 1e10 / A[nr, nphi]
energymin = 1e-10*energy0
energy = np.ones((Nrad, Naz), dtype=np.float64)*energymin

energy[nr, nphi] = energy0

# save this energy array to file
energy.tofile("output/out/snapshots/0/energy.dat")

import os
if os.path.exists("output/out/snapshots/damping"):
    energy.tofile("output/out/snapshots/damping/energy.dat")

