#!/usr/bin/env python3

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const


def main():
    fig = plot_drift("output")
    fig.savefig("trajectories.jpg")


def plot_drift(outdir):

    fig, ax = plt.subplots()

    particles = construct_dust_trajectories(outdir)

    sizes = []
    Sts = []
    rdot_abs = []

    for i, p in particles.items():
        # if not i%10 == 0:
        #     continue

        t = p["time"]
        r = p["r"]
        
        
        size = p["size"][0].to("cm")
        stokes = p["stokes"]
        vtheo = vdrift_theo(stokes, r.to("au"))
        print(stokes)
        vtheo = vtheo.to("cm/s")
        
        r = r.to("cm")
        t = t.to("s")
        rdot = (r[1:] - r[:-1])/(t[1:] - t[:-1])
        rdot = rdot.to("cm/s")
        
        print(rdot[-1]/vtheo[-1])
        
        line, = ax.plot(t[:-1].to("yr"), -rdot, label=f"s = {size:.1e}")
        color = line.get_color()
        ax.plot(t.to("yr"), -vtheo, ls=":", color=color)

    ax.set_yscale("log")
    ax.set_ylabel(r"$-\dot{r}$ [cm/s]")
    ax.set_xlabel(r"$t$ [orbits]")
    
    ax.legend()

    print(vdrift_theo(7.4e-1, 1*u.au))
    print(vdrift_theo(7.4e-1, 1*u.au).to("au/yr"))


    return fig


def vdrift_theo(stokes, r):
    """ Drift speed according to Picogna & Kley 2015 Eq. (C.1) (10.1051/0004-6361/201526921)
    and Nakagawa+1986 Eq. (1.9) (10.1016/0019-1035(86)90121-1). 
    Note that in Eq. (1.9) for eta, there is a '/' missing and r OmegaK^2 needs to be in the denominator.  
    The implementation matches the values in the plot of Picogna & Kley 2015 Fig. C.2 
    (though the mislabeled the unit on the y axis to be cm/s but it must be au/yr)."""
    Mstar = 1*u.solMass
    h = 0.05
    # r = 1*u.au
    vK = np.sqrt(const.G*Mstar/r).decompose()
    eta = 0.5*h**2
    vdrift = -eta*vK/(stokes + stokes**-1)
    return vdrift
    
    

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
