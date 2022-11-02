#!/usr/bin/env python3

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from drift_theo import vdrift_theo

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
            
            # print("stokes", stokes[-1], "\t size = ", size, "\t sim / theo =", rdot[-1]/vtheo[-1])
            
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
    # print("time", t)
    # print("timeunit", timeunit)
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
        phi = particles["phi"].value*u.rad

        x = r*np.cos(phi)
        y = r*np.sin(phi)

    particles["x"] = x
    particles["y"] = y
    particles["r"] = r
    particles["phi"] = phi

    return particles


if __name__ == "__main__":
    main()
