#!/usr/bin/env python3

import argparse

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from simdata import Data


def parse_cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir", help="Path to the output dir.")
    parser.add_argument("Ns", type=int, nargs="+", help="Output number")
    parser.add_argument("--extradirs", nargs="+", type=str, default=[], help="Paths to extra output dirs.")
    parser.add_argument("--outfile", help="Filename of the output file.")
    parser.add_argument("--labels", type=str, nargs="+", help="Labels identifying the different output dirs.")
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

        data = Data(outDir)

        for N in Ns:
            try:
                plot(ax, data, N, name=label)
            except KeyError:
                print(f"Can't plot output step {N}")

    if args.outfile:
        fig.savefig(args.outfile, dpi=300)
    else:
        plt.show()


def plot(ax, data, N, name=""):

    print(f"Plotting particle distribution for N = {N}")

    sizes = [1e-5*u.cm]

    for size in sizes:
        particles, time = data.particles.get(N)
        r = particles["r"].to_value("au")
        N_particles = len(r)
        # hacky way around fargocpt outputting x and y
        # for adaptive integrator
        if any(r < 0):
            phi = particles["phi"].value
            r = np.sqrt(r**2 + phi**2)
        print("N = {} particles".format(len(r)))
        label = "{} t = {:.2f}, N = {}".format(name, time.to("yr"), N_particles)
        counts, interfaces = np.histogram(r, bins=51)
        mid = 0.5*(interfaces[1:] + interfaces[:-1])
        dr = interfaces[1:] - interfaces[:-1]
        # sigma_dust = counts/dr*mid
        # sigma_dust /= np.max(sigma_dust)
        sigma_dust = counts*1
        # sigma_dust = sigma_dust / np.max(sigma_dust)
        ax.plot(mid, sigma_dust, label=label)

    ax.set_xlabel(r"$r$ [au]")
    ax.set_ylabel(r"particle histogram normalized")

    # ax.set_xscale("log")
    ax.set_yscale("log")

    ax.legend(loc="lower right")

    title = "Particle location histogram"
    ax.set_title(title)

    ax.grid(True, alpha=0.6)
    # ax.grid(axis="x", which="minor", alpha=0.6)

    ax.axvline(13.54, alpha=0.75, color="black", ls=":")
    ax.axvline(7.17, alpha=0.75, color="black", ls=":")

    ax.set_ylim(bottom=1, top=1000)

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


if __name__ == "__main__":
    main()
