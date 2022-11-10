#!/usr/bin/env python3

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from scipy import interpolate

from drift_theo import vdrift_theo
from load_dust import construct_dust_trajectories

def main():
    fig = plot_drift("output/dust_drift")
    fig.savefig("drift.jpg")


def plot_drift(outdir):

    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(3, 1)

    ax = fig.add_subplot(gs[:2, :])

    particles = construct_dust_trajectories(outdir)

    sizes = []
    Sts = []
    rdot_abs = []

    for i, p in particles.items():

        t = p["time"].to_value("s")
        r = p["r"].to_value("cm")

        rdot = np.abs((r[-1] - r[-2])/(t[-1] - t[-2]))
        # rdot = (rdot*u.cm/u.s).to_value(5.2*u.au/(5.2**1.5*u.yr))
        size = p["size"][0].to_value("cm")
        St = p["stokes"][0]

        rdot_abs.append(rdot)
        sizes.append(size)
        Sts.append(St)

    rdot_abs = np.array(rdot_abs)
    sizes = np.array(sizes)
    Sts = np.array(Sts)

    # X = sizes
    X = Sts
    Y = rdot_abs
    ax.plot(X, Y, marker="x")

    ax.set_ylabel(r"$-\dot{r}$ [cm/s]")

    ax.set_xscale("log")
    ax.set_yscale("log")

    # secondary y axis
    secax = ax.secondary_yaxis(
        'right', functions=(cmps_to_aupyr, aupyr_to_cmps))
    secax.set_ylabel(r"$\dot{r}$ [au/yr]")

    ax.grid(alpha=0.3)
    
    vdrift = vdrift_theo(Sts, 1*u.au).to_value("cm/s")
    ax.plot(Sts, -vdrift)

    try:
        f = interpolate.InterpolatedUnivariateSpline(sizes, Sts)
        f2 = interpolate.InterpolatedUnivariateSpline(Sts, sizes)
        secax = ax.secondary_xaxis(
            'top', functions = (f2, f)
        )
        secax.set_xlabel("size [cm]")
    except Exception:
        pass


    ax = fig.add_subplot(gs[2:3, :])
    reldiff = (np.abs(rdot_abs - np.abs(vdrift))) / np.abs(vdrift)
    ax.plot(Sts, reldiff)
    ax.set_xscale("log")
    ax.set_ylabel("rel. diff.")
    ax.set_xlabel(r"Stokes number")
    ax.set_ylim(bottom=0)
    ax.grid()


    # ax = fig.add_subplot(gs[3:])
    # ax.plot(sizes, Sts)
    # ax.set_yscale("log")
    # ax.set_xscale("log")
    # ax.grid()
    # ax.set_xlabel(r"$s$ [cm]")
    # ax.set_ylabel("Stokes number")

    return fig
    

def cmps_to_aupyr(x):
    factor = (1*u.au/u.yr).to("cm/s").value
    return x / factor


def aupyr_to_cmps(x):
    factor = (1*u.au/u.yr).to("cm/s").value
    return x * factor


if __name__ == "__main__":
    main()
