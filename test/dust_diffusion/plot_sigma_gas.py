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
    args = parser.parse_args()
    return args


def main():
    args = parse_cli_args()
    Ns = args.Ns
    outDir = args.outdir

    data = Data(outDir)

    fig, ax = plt.subplots()
    for N in Ns:
        plot(ax, data, N)

    plot_name = "particle_histogram"
    plt.show()


def plot(ax, data, N):

    print(f"Plotting 1d surface density for N = {N}")

    sig = data.fluids["gas"].get("1d", "mass density", N)
    X = sig.grid.get_coordinates("r").to_value("au")
    Y = sig.data.to_value("g/cm2")

    label = "t = {:.2f}".format(sig.time.to("yr"))

    ax.plot(X, Y, label=label)

    ax.set_xlabel(r"$r$ [au]")
    ax.set_ylabel(r"$\Sigma$ [g/cm2]")

    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # dx = 1e-5
    # ax.set_xlim(left=10-dx, right=10+dx)
    # ax.set_ylim(bottom=0.5, top=5500)

    ax.grid(alpha=0.3)
    # ax.grid(axis="x", which="minor", alpha=0.3)

    ax.legend(loc="upper left")

    title = "Surface density"
    ax.set_title(title)


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
