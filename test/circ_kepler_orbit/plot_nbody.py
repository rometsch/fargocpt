#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as units
from fargocpt import Loader

def main():
    l = Loader("../../output/tests/circ_kepler_orbit/out")

    fig, axes = plt.subplots(nrows=2, ncols=2)

    Npl = len(l.nbody)
    for n in range(Npl):
        nb = l.nbody[n]
        a = nb.semi_major_axis.to_value("au")
        t = nb.time.to_value("kyr")
        ax = axes[0,0]
        ax.plot(t, np.abs(a-1), label=f"{n}")

        e = nb.eccentricity.value
        ax = axes[0,1]
        ax.plot(t, e)

    ax = axes[0,0]
    ax.legend()
    ax.set_yscale("log")
    ax.set_ylabel(r"|$a$-1au| / 1au")

    ax = axes[0,1]
    ax.set_ylabel(r"$e$")
    ax.set_xlabel(r"$t$ [kyr]")

    for n in range(Npl):
        nb = l.nbody[n]
        t = nb.time.to_value("kyr")
        x = nb.x.to_value("au")
        ax = axes[1,0]
        ax.plot(t, x, label=f"{n}")

        nb = l.nbody[n]
        t = nb.time.to_value("kyr")
        y = nb.y.to_value("au")
        ax = axes[1,1]
        ax.plot(t, y)

    ax = axes[1,0]
    ax.legend()
    ax.set_ylabel(r"$x$ [au]")

    ax = axes[1,1]
    ax.set_ylabel(r"$y$ [au]")
    ax.set_xlabel(r"$t$ [kyr]")



    fig.savefig("plot.jpg")

    plt.show()

if __name__=="__main__":
    main()