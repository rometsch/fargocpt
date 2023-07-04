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

ri = np.linspace(rmin, rmax, Nrad+1)
phii = np.linspace(0, 2*np.pi, Naz+1)
Ri, Phii = np.meshgrid(ri, phii, indexing="ij")
rc = 2/3*(ri[1:]**2/(ri[1:]+ri[:-1])) + ri[:-1] # approx center in polar coords
phic = 0.5*(phii[1:]+phii[:-1])

Rc, Phic = np.meshgrid(rc, phic, indexing="ij")
Xc = Rc*np.cos(Phic)
Yc = Rc*np.sin(Phic)

dphi = phii[1] - phii[0]
dr = ri[1:] - ri[:-1]
A = 0.5*(Ri[1:,1:]**2 - Ri[:-1,1:]**2)*dphi

# now create an energy array with values in cgs with the energy located in one cell

nr = Nrad//2
nphi = 0

# print(dr[nr])

c = 1e10
rho = 1
kappa = 1
K = 1/3*c/(rho*kappa)

xcell = Xc[nr, nphi]
ycell = Yc[nr, nphi]
DX = Xc - xcell
DY = Yc - ycell
Dist = np.sqrt(DX**2 + DY**2)

def analytical_solution(r, t, E0, K):
    return 3*E0/(4*np.pi*t*K)*np.exp(-r**2/(4*K*t))

t0 = 1e-10

energy0 = 1e10 / A[nr, nphi]
# energymin = 1e-10*energy0
# energymin = 0
# energy = np.ones((Nrad, Naz), dtype=np.float64)*energymin

initial_condition = analytical_solution(Dist,t0, energy0, K)
energy = initial_condition

# save this energy array to file
energy.tofile("output/out/snapshots/0/energy.dat")

import os
if os.path.exists("output/out/snapshots/damping"):
    energy.tofile("output/out/snapshots/damping/energy.dat")

