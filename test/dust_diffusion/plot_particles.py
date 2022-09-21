#!/usr/bin/env python3

import argparse

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from disgrid import Data


def parse_cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir", help="Path to the output dir.")
    parser.add_argument("N", type=int, help="Output number")
    parser.add_argument("--outfile", help="Filename of the output file.")
    args = parser.parse_args()
    return args


def main():
    args = parse_cli_args()
    N = args.N

    outDir = args.outdir
    
    fig, ax = plt.subplots()
    
    data = Data(outDir)

    plot(ax, data, N)
            
    if args.outfile:
        fig.savefig(args.outfile, dpi=300)
    else:
        plt.show()


def plot(ax, data, N, name=""):

    print(f"Plotting particle distribution for N = {N}")
    particles, time = data.particles.get(N)

    # hacky way around fargocpt outputting x and y
    # for adaptive integrator
    cartesian = any(particles["r"] < 0)
    if cartesian:
        x = particles["r"]
        y = particles["phi"]
    else:
        r = particles["r"]
        phi = particles["phi"]
    
        x = r*np.cos(phi)
        y = r*np.sin(phi)



    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    title = f"Particle distribution N={N}"
    ax.set_title(title)

    ax.grid(True, alpha=0.6)
    
    ax.scatter(x, y)
    ax.set_aspect("equal")
    
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
