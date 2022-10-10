#!/usr/bin/env python3

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


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

    plot(ax, outDir, N)
            
    if args.outfile:
        fig.savefig(args.outfile, dpi=300)
    else:
        plt.show()


def plot(ax, datadir, N):

    print(f"Plotting particle distribution for N = {N}")
    particles = get_particles(datadir, N)

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
