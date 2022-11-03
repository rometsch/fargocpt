#!/usr/bin/env python3

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from drift_theo import vdrift_theo
from load_dust import construct_dust_trajectories

def main():
    fig = plot_trajectory("output/dust_drift")
    fig.savefig("trajectories.jpg")


def plot_trajectory(outdir):

    fig, axes = plt.subplots(figsize=(10,5), ncols=2)

    particles = construct_dust_trajectories(outdir)

    sizes = []
    Sts = []
    rdot_abs = []

    for i, p in particles.items():
        try:
            # if not i%10 == 0:
            #     continue

            t = p["time"]
            r = p["r"]
            
            
            size = p["size"][0].to("cm")
            stokes = p["stokes"]
            vtheo = vdrift_theo(stokes, r.to("au"))
            vtheo = vtheo.to("cm/s")
            
            # print(t, r, stokes)
            
            r = r.to("cm")
            t = t.to("s")
            rdot = (r[1:] - r[:-1])/(t[1:] - t[:-1])
            rdot = rdot.to("cm/s")
            
            print("stokes", stokes[-1], "\t size = ", size, "\t sim / theo =", rdot[-1]/vtheo[-1])
            
            ax = axes[0]
            line, = axes[0].plot(t[:-1].to("yr"), -rdot, label=f"s = {size:.1e}")
            color = line.get_color()
            ax.plot(t.to("yr"), -vtheo, ls=":", color=color)
            
            axes[1].plot(t.to("yr"), r.to("au"), color=color, label=f"s = {size:.1e}")
        except IndexError as e:
            print(e)

    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$-\dot{r}$ [cm/s]")
    axes[0].set_xlabel(r"$t$ [orbits]")
    
    axes[1].legend()

    # print(vdrift_theo(7.4e-1, 1*u.au))
    # print(vdrift_theo(7.4e-1, 1*u.au).to("au/yr"))


    return fig


def cmps_to_aupyr(x):
    factor = (1*u.au/u.yr).to("cm/s").value
    return x / factor


def aupyr_to_cmps(x):
    factor = (1*u.au/u.yr).to("cm/s").value
    return x * factor



if __name__ == "__main__":
    main()
