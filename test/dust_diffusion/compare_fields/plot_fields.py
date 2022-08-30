#!/usr/bin/env python3

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from simdata import Data


def parse_cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fieldname", type=str, help="Name of the variable to plot.")
    parser.add_argument("outdir", help="Path to the output dir.")
    parser.add_argument("Ns", type=int, nargs="+", help="Output number")
    parser.add_argument("--extradirs", nargs="+", type=str, default=[], help="Paths to extra output dirs.")
    parser.add_argument("--outfile", help="Filename of the output file.")
    args = parser.parse_args()
    return args


def main():
    args = parse_cli_args()
    Ns = args.Ns

    outDirs = [args.outdir] + args.extradirs
    
    fig, ax = plt.subplots()
    
    for outDir in outDirs:

        for N in Ns:
            try:
                plot(ax, outDir, args.fieldname, N)
            except KeyError:
                print(f"Can't plot output step {N}")

    if args.outfile:
        fig.savefig(args.outfile, dpi=300)
    else:
        plt.show()


def plot(ax, outDir, varname, N):

    datafile = os.path.join(outDir, f"gas{varname}{N}.dat")

    rad_file = os.path.join(outDir, "used_rad.dat")
    rs = np.genfromtxt(rad_file)
    N_rad = len(rs)-1
    rs = 0.5*(rs[1:] + rs[:-1])

    data = np.fromfile(datafile)
    N_cells = len(data)
    N_az = int(N_cells/N_rad)

    data_2d = data.reshape(N_rad, N_az)

    data_1d_avg = np.average(data_2d, axis=1)

    label = outDir
    ax.plot(rs, data_1d_avg, label=label)

    ax.set_xlabel(r"$r$ [au]")
    ax.set_ylabel(r"1d average [code units]")

    # ax.set_xscale("log")
    # ax.set_yscale("log")

    ax.legend(loc="lower right")

    title = varname
    ax.set_title(title)

    ax.grid(True, alpha=0.6)
    # ax.grid(axis="x", which="minor", alpha=0.6)

    ax.axvline(13.54, alpha=0.75, color="black", ls=":")
    ax.axvline(7.17, alpha=0.75, color="black", ls=":")

    # ax.set_ylim(bottom=1, top=1000)


if __name__ == "__main__":
    main()
