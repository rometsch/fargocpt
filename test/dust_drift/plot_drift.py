#!/usr/bin/env python3

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from scipy import interpolate

from drift_theo import vdrift_theo

from matplotlib.ticker import NullLocator

from fargocpt import Loader

def main():
    fig = plot_drift("../../output/tests/dust_drift/out")
    fig.savefig("drift.jpg", dpi=150, bbox_inches="tight")
    # fig.savefig("drift.pdf", dpi=300, bbox_inches="tight")


def plot_drift(outdir):


    fig, axs = plt.subplots(2,1,height_ratios=[3,1], sharex=True, figsize=(6,4))
    fig.subplots_adjust(hspace=0)

    ax = axs[0]
    l = Loader(outdir)
    particles = l.particles.timeseries(["r", "stokes", "size"])

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

    vdrift = vdrift_theo(Sts, 1*u.au).to_value("cm/s")


    # X = sizes
    X = Sts
    Y = rdot_abs
        
    ax.plot(X, Y, lw=2, label="code", marker="+", markersize=10, ls="-")
    ax.plot(Sts, -vdrift, lw=2, ls="--", label="theoretical")
    
    

    ax.set_ylabel(r"$-\dot{r}$ [cm/s]")

    ax.set_xscale("log")
    ax.set_yscale("log")
    
    # secondary y axis
    secax = ax.secondary_yaxis(
        'right', functions=(cmps_to_aupyr, aupyr_to_cmps))
    secax.set_ylabel(r"$\dot{r}$ [au/yr]")
    secax.set_yticks([1e-8, 1e-6, 1e-4, 1e-2])
    secax.yaxis.set_minor_locator(NullLocator())

    # ax.grid(alpha=0.3)
    
    ax.set_yticks([1e-2, 1e0, 1e2, 1e4])
    ax.yaxis.set_minor_locator(NullLocator())



    ax.legend()

    # try:
    #     f = interpolate.InterpolatedUnivariateSpline(sizes, Sts)
    #     f2 = interpolate.InterpolatedUnivariateSpline(Sts, sizes)
    #     secax = ax.secondary_xaxis(
    #         'top', functions = (f2, f)
    #     )
    #     secax.set_xlabel("size [cm]")
    # except Exception:
    #     pass

    ax = axs[1]
    reldiff = (np.abs(rdot_abs - np.abs(vdrift))) / np.abs(vdrift)
    ax.plot(Sts, reldiff, lw=2)
    ax.set_xscale("log")
    ax.set_ylabel("rel. diff.")
    ax.set_xlabel(r"Stokes number")
    ax.set_ylim(bottom=1e-4, top=2)
    ax.set_yscale("log")
    # ax.grid()


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
