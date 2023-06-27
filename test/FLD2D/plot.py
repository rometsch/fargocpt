#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("N", type=int, default=-1, help="Snapshot number to plot.")
parser.add_argument("-n", type=str, default="energy density", help="Variable to plot")
parser.add_argument("-p", action="store_true", help="Print the values")
parser.add_argument("--log", action="store_true", help="Log plot")


opts = parser.parse_args()

from disgrid import Data
d = Data("output/out")
d.avail()


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors


standard_names = [
    "mass density",
    "energy density",
    "velocity radial",
    "velocity azimuthal"
]

def plot_field(data, name, N, ax=None, dataunit=None, lengthunit=None, vmin=None, vmax=None, cmap="viridis", diff=False, toprint=False, logplot=False):
    
    if name in standard_names:
        field = data.get(var=name, dim="2d", N=N)
        if dataunit is None:
            dataunit = field.data.unit
        Z = field.data.to_value(dataunit)
        if diff:
            Z0 = data.get(var=name, dim="2d", N=0).data.to_value(dataunit)
            Z = Z - Z0
    else:
        field = data.get(var="mass density", dim="2d", N=N)
        if dataunit == None:
            dataunit = "a.u."

    ri = field.grid.get_interfaces("r")
    if lengthunit is None:
        lengthunit = ri.unit
    ri = ri.to_value(lengthunit)
    phii = field.grid.get_interfaces("phi").to_value("rad")
    PHI, R = np.meshgrid(phii, ri)
    X = R*np.cos(PHI)
    Y = R*np.sin(PHI)

    Nr = len(ri)-1
    Nphi = len(phii)-1

    if name not in standard_names:
        
        Z = np.fromfile(data.path+f"/snapshots/{N}/{name}.dat", dtype=np.float64)
        # handle quantities defined on interfaces
        try:
            Z = Z.reshape(Nr, Nphi)
        except ValueError:
            Z = Z.reshape(Nr+1, Nphi)[:-1,:]
        if diff:
            Z0 = np.fromfile(data.path+f"/snapshots/{0}/{name}.dat", dtype=np.float64)
            try:
                Z0 = Z0.reshape(Nr, Nphi)
            except ValueError:
                Z0 = Z0.reshape(Nr+1, Nphi)[:-1,:]
            Z = Z-Z0

    if not diff:
        print("min value: {:e}".format(np.min(Z)))
        print("max value: {:e}".format(np.max(Z)))


    if toprint:
        print(Z)

    if ax is None:
        fig, ax = plt.subplots(dpi=150)
    else:
        fig = ax.get_figure()

    if logplot:
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
    if diff:
        vmax = np.max(np.abs(Z))
        norm = mplcolors.Normalize(vmin=-vmax, vmax=vmax)
        cmap = "bwr"
    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    t = field.time.to_value("s")
    ax.set_title(f" t={t:.2e}s, N={N}")

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(f"{name} [{dataunit}]")

    ax.set_xlabel(f"x [{lengthunit}]")
    ax.set_ylabel(f"y [{lengthunit}]")
    
    return fig

if opts.N < 0:
    N = d.avail()["Nlast"] + opts.N + 1
else:
    N = opts.N

varname = opts.n

try:
    E0 = d.get(var=varname, dim="2d", N=0).data.to_value("erg/cm2")
    Elast = d.get(var=varname, dim="2d", N=N).data.to_value("erg/cm2")
    print("Difference to initial =", np.sum(np.abs(E0 - Elast)))
except KeyError:
    pass


plot_field(d, varname, N, dataunit=None, cmap="magma", diff=False, toprint=opts.p, logplot=opts.log);
plot_field(d, varname, N, dataunit=None, cmap="magma", diff=True, toprint=opts.p);
plt.show()
