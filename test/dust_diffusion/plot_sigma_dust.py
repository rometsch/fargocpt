#!/usr/bin/env python3

import argparse

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

import os

def parse_cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir", help="Path to the output dir.")
    parser.add_argument("Ns", type=int, nargs="+", help="Output number")
    parser.add_argument("--extradirs", nargs="+", type=str, default=[], help="Paths to extra output dirs.")
    parser.add_argument("--outfile", help="Filename of the output file.")
    parser.add_argument("--labels", type=str, nargs="+", help="Labels identifying the different output dirs.")
    parser.add_argument("--normalize", action="store_true", help="Normalize the distribution to the maximum value.")
    args = parser.parse_args()
    return args


def main():
    args = parse_cli_args()
    Ns = args.Ns

    outDirs = [args.outdir] + args.extradirs
    if args.labels is not None:
        labels = args.labels
    else:
        labels = [""] * len(outDirs)
    
    fig, ax = plt.subplots()
    
    for label, outDir in zip(labels, outDirs):

        for N in Ns:
            try:
                plot(ax, outDir, N, name=label, normalize=args.normalize)
            except KeyError:
                print(f"Can't plot output step {N}")

    if args.outfile:
        fig.savefig(args.outfile, dpi=300)
    else:
        plt.show()


def get_closest_analytic(t):
    profiles = np.load("analytic_profiles.npy")
    times = np.load("analytic_times.npy")
    r = np.load("analytic_radii.npy")
    
    n = np.argmin(np.abs(times - t))
    return r, profiles[n], times[n]

def plot(ax, outdir, N, name="", normalize=False):

    print(f"Plotting particle distribution for N = {N}")

    particles = get_particles(outdir, N)
    r = particles["r"]
    N_particles = len(r)
    # hacky way around fargocpt outputting x and y
    # for adaptive integrator
    if any(r < 0):
        phi = particles["phi"]
        r = np.sqrt(r**2 + phi**2)
    print("N = {} particles".format(len(r)))
    t = get_time(outdir, N)
    
    label = "{} N = {}, t = {:.5g} yr, Nparts = {}".format(name, N, t, N_particles)
    counts, interfaces = np.histogram(r, bins=51)
    mid = 0.5*(interfaces[1:] + interfaces[:-1])
    dr = interfaces[1:] - interfaces[:-1]
    sigma_dust = counts/(dr*mid)
    # sigma_dust /= np.max(sigma_dust)
    # sigma_dust = counts*1
    if normalize:
        sigma_dust = sigma_dust / np.max(sigma_dust)
    line, = ax.plot(mid, sigma_dust, label=label)

    ra, pa, ta = get_closest_analytic(t)
    pa = pa/np.max(pa) *np.max(sigma_dust)
    ax.plot(ra, pa, color=line.get_color(), ls="--")#, label=f"analytic t = {ta:3g} yr")

    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"particle histogram normalized")

    # ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend(loc="best")

    title = "Particle location histogram"
    ax.set_title(title)

    ax.grid(True, alpha=0.6)
    # ax.grid(axis="x", which="minor", alpha=0.6)

    # ax.axvline(13.54, alpha=0.75, color="black", ls=":")
    # ax.axvline(7.17, alpha=0.75, color="black", ls=":")

    # ax.set_ylim(bottom=1e-5, top=1)
    if normalize:
        ax.set_ylim(bottom=1e-3, top=1)
    else:
        ax.set_ylim(bottom=1, top=1e5)
        
    ax.set_xlim(left=1, right=20)

def particles_by_size(data, N, size=[-np.inf, np.inf]):
    """ Construct trajectories for all dust particles. """
    try:
        len(size)
    except TypeError:
        size = [0.9*size, 1.1*size]

    vals, time = data.particles.get(N)
    ids = vals["id"]
    sizes = vals["size"]

    mask = np.logical_and(sizes >= size[0], sizes <= size[1])
    for key in vals:
        vals[key] = vals[key][mask]

    return vals, time

def get_time(datadir, N):
    filename = os.path.join(datadir, f"snapshots/timeSnapshot.dat")
    times = np.genfromtxt(filename, usecols=(2))
    units = get_units(datadir)
    times = (times*units["time"]).to_value("yr")
    return times[N]

def get_units(datadir):
    filename = os.path.join(datadir, "units.dat")
    units = {}
    with open(filename, "r") as infile:
        for line in infile:
            line = line.strip()
            if line == "" or line[0] == "#":
                continue
            else:
                parts = line.split()
                if len(parts) == 3:
                    units[parts[0]] = parts[1] * u.Unit(parts[2])
    return units

def get_particles(datadir, N):
    filename = os.path.join(datadir, f"snapshots/{N}/particles.dat")
    res = np.fromfile(
        filename, dtype=[('Id', np.dtype(int)), ('Values', np.dtype(float), 11)])
    ids = res["Id"]
    vals = res["Values"]

    particles = {
        "id": ids,
        "r": vals[:, 0],
        "phi": vals[:, 1],
        "r dot": vals[:, 2],
        "phi dot": vals[:, 3],
        "r ddot": vals[:, 4],
        "phi ddot": vals[:, 5],
        "mass": vals[:, 6],
        "size": vals[:, 7],
        "timestep": vals[:, 8],
        "facold": vals[:, 9],
        "stokes": vals[:, 10]
    }
    return particles


if __name__ == "__main__":
    main()
