#!/usr/bin/env python3

import numpy as np
import yaml


# get size of grid specified in the setup file
with open("setup.yml", "r") as infile:
    params = yaml.safe_load(infile)
    Nrad = params["Nrad"]
    Naz = params["Nsec"]


# now create an energy array with values in cgs with the energy located in one cell
energymin = 1
energy = np.ones((Nrad, Naz), dtype=np.float64)*energymin
energy0 = 1e8

Ncells = 5
# energy[:,:] = energy0/2
energy[3*Nrad//4:3*Nrad//4+Ncells, 0:Ncells] = energy0

# save this energy array to file
energy.tofile("output/out/snapshots/0/energy.dat")

import os
if os.path.exists("output/out/snapshots/damping"):
    energy.tofile("output/out/snapshots/damping/energy.dat")

