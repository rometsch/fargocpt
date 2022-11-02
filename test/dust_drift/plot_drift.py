#!/usr/bin/env python3

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from scipy import interpolate

from drift_theo import vdrift_theo

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

    f = interpolate.InterpolatedUnivariateSpline(sizes, Sts)
    f2 = interpolate.InterpolatedUnivariateSpline(Sts, sizes)
    secax = ax.secondary_xaxis(
        'top', functions = (f2, f)
    )
    secax.set_xlabel("size [cm]")



    ax = fig.add_subplot(gs[2:3, :])
    reldiff = (rdot_abs - np.abs(vdrift)) / np.abs(vdrift)
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


def get_output_times(datadir):
    filepath = os.path.join(datadir, "snapshots/timeSnapshot.dat")
    t = np.genfromtxt(filepath, usecols=2)
    # get unit
    timeunit = None
    with open(filepath, "r") as infile:
        for line in infile:
            m = re.search("physical time \| ([0-9e+-\.]+) s", line)
            if m is not None:
                timeunit = u.Unit(m.groups()[0] + "s")
    if timeunit is not None:
        t = t*timeunit
    return t


def construct_dust_trajectories(datadir):
    """ Construct trajectories for all dust particles. """
    particles = {}
    times = get_output_times(datadir)

    for n, time in enumerate(times):
        vals = get_particles(datadir, n)
        ids = vals["id"]

        if n == 0:
            for i in ids:
                particles[i] = {}
                for key in [k for k in vals.keys() if not k == "id"]:
                    particles[i][key] = []
                particles[i]["time"] = []

        for key in [k for k in vals.keys() if not k == "id"]:
            for i, v in zip(ids, vals[key]):
                particles[i][key].append(v)
        for i in ids:
            particles[i]["time"].append(time)

    for ind, p in particles.items():
        for varname in p:
            p[varname] = u.Quantity(p[varname])

    return particles


def get_particles(datadir, N):
    filename = os.path.join(datadir, f"snapshots/{N}/particles.dat")
    res = np.fromfile(
        filename, dtype=[('Id', np.dtype(int)), ('Values', np.dtype(float), 11)])
    ids = res["Id"]
    vals = res["Values"]

    L0 = u.Unit("1 au")

    particles = {
        "id": ids,
        "r": vals[:, 0]*L0,
        "phi": vals[:, 1]*L0,
        "r dot": vals[:, 2],
        "phi dot": vals[:, 3],
        "r ddot": vals[:, 4],
        "phi ddot": vals[:, 5],
        "mass": vals[:, 6],
        "size": vals[:, 7]*L0,
        "timestep": vals[:, 8],
        "facold": vals[:, 9],
        "stokes": vals[:, 10]
    }

    # hacky way around fargocpt outputting x and y
    # for adaptive integrator
    is_cartesian = any(particles["r"] < 0)
    if is_cartesian:
        x = particles["r"]
        y = particles["phi"]

        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(x, y)
    else:
        r = particles["r"]
        phi = particles["phi"].value

        x = r*np.cos(phi)
        y = r*np.sin(phi)

    particles["x"] = x
    particles["y"] = y
    particles["r"] = r
    particles["phi"] = phi

    return particles


if __name__ == "__main__":
    main()
