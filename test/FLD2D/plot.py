#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("N", type=int, default=-1, help="Snapshot number to plot.")
parser.add_argument("-n", type=str, default="energy density", help="Variable to plot")
parser.add_argument("-p", action="store_true", help="Print the values")
parser.add_argument("--log", action="store_true", help="Log plot")
parser.add_argument("-d", type=str, default=None, help="Output directory")
parser.add_argument("--diff", action="store_true", help="Print difference to first.")


opts = parser.parse_args()

outdir = opts.d if opts.d is not None else "output/out"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors

from create_input import get_fargo_grid

def plot_field(outdir, name, N, ax=None, dataunit=None, lengthunit=None, vmin=None, vmax=None, cmap="viridis", diff=False, toprint=False, logplot=False):
    
    setupfile = outdir + f"/snapshots/{N}/config.yml"
    g = get_fargo_grid(setupfile)

    X = g.Xi
    Y = g.Yi

    Nr = len(g.ri)-1
    Nphi = len(g.phii)-1
        
    Z = np.fromfile(outdir+f"/snapshots/{N}/{name}.dat", dtype=np.float64)
    # handle quantities defined on interfaces
    try:
        Z = Z.reshape(Nr, Nphi)
    except ValueError:
        Z = Z.reshape(Nr+1, Nphi)[:-1,:]
    if diff:
        Z0 = np.fromfile(outdir+f"/snapshots/{0}/{name}.dat", dtype=np.float64)
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
        if (Z<=0).all():
            Z = -Z
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
    if diff:
        vmax = np.max(np.abs(Z))
        if logplot:
            norm = mplcolors.SymLogNorm(vmin=-vmax, vmax=vmax, linthresh = 1e-2*vmax)
        else:
            norm = mplcolors.Normalize(vmin=-vmax, vmax=vmax)
        cmap = "bwr"
    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    ax.set_title(f"N={N}")

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(f"{name} [{dataunit}]")

    ax.set_xlabel(f"x [{lengthunit}]")
    ax.set_ylabel(f"y [{lengthunit}]")
    
    return fig

N = opts.N

varname = opts.n



fig = plot_field(outdir, varname, N, dataunit=None, cmap="magma", diff=False, toprint=opts.p, logplot=opts.log)
fig.savefig("plot.jpg", dpi=150)
if opts.diff:
    plot_field(outdir, varname, N, dataunit=None, cmap="magma", diff=True, toprint=opts.p, logplot=opts.log)
plt.show()
