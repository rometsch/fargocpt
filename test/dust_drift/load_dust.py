#!/usr/bin/env python3

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

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

    # import yaml
    # outdir = 'output/out'
    # with open(outdir + '/infoParticles.yml', 'r') as infile:
    #     infoParticles = yaml.safe_load(infile)

    # N = 0
    # filename = f'output/out/snapshots/{N}/particles.dat'

    # keys = [k for k in infoParticles['variables'].keys()]
    # types = []
    # for key, info in infoParticles['variables'].items():
    #     dtype = info['type']
    #     if dtype == 'double':
    #         types.append((key, np.dtype(np.dtype('f8'))))
    #     elif dtype == 'unsigned long':
    #         types.append((key, 'u8'))
    #     else:
    #         raise ValueError(f'Unknown type {dtype}')

    # res = np.fromfile(filename, dtype=types)



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