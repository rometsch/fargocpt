#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("N", type=int, default=-1, help="Snapshot number to plot.")
parser.add_argument("dir1", type=str, default=None, help="Output directory nr 1")
parser.add_argument("dir2", type=str, default=None, help="Output directory nr 1")
parser.add_argument("-n", "--name", type=str, default="energy density", help="Variable to plot")
parser.add_argument("-p", action="store_true", help="Print the values")
parser.add_argument("--log", action="store_true", help="Log plot")


opts = parser.parse_args()

from disgrid import Data

datas = [Data(opts.dir1),
         Data(opts.dir2)]


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors

N = opts.N
logplot = opts.log

field = datas[0].get(var="mass density", dim="2d", N=N)

ri = field.grid.get_interfaces("r")
lengthunit = ri.unit

ri = ri.to_value(lengthunit)
phii = field.grid.get_interfaces("phi").to_value("rad")
PHI, R = np.meshgrid(phii, ri)
X = R*np.cos(PHI)
Y = R*np.sin(PHI)

Nr = len(ri)-1
Nphi = len(phii)-1

name = opts.name
Z1 = np.fromfile(datas[0].path+f"/snapshots/{N}/{name}.dat", dtype=np.float64).reshape(Nr, Nphi)
Z2 = np.fromfile(datas[1].path+f"/snapshots/{N}/{name}.dat", dtype=np.float64).reshape(Nr, Nphi)
Z = Z2-Z1
Z = Z / np.max(Z1)

print("min value: {:e}".format(np.min(Z)))
print("max value: {:e}".format(np.max(Z)))

fig, ax = plt.subplots(dpi=150)

vmax = np.max(np.abs(Z))
if logplot:
    norm = mplcolors.SymLogNorm(vmin=-vmax, vmax=vmax, linthresh = 1e-2*vmax)
else:
    norm = mplcolors.Normalize(vmin=-vmax, vmax=vmax)
cmap = "bwr"
pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
ax.set_aspect("equal")

t = field.time.to_value("s")
ax.set_title(f" t={t:.2e}s, N={N}")

cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label(f"{name}")

ax.set_xlabel(f"x [{lengthunit}]")
ax.set_ylabel(f"y [{lengthunit}]")

plt.show()