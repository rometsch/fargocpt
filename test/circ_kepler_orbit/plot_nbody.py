#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as units
import disgrid

def main():
    d = disgrid.Data("../../output/tests/circ_kepler_orbit/out")

    fig, axes = plt.subplots(nrows=2, ncols=2)

    Npl = len(d.avail()["planets"])
    for n in range(Npl):
        q = d.get("semi-major axis", planet=n)
        t = q.time.to_value("kyr")
        a = q.data.to_value("au")
        ax = axes[0,0]
        ax.plot(t, np.abs(a-1), label=f"{n}")

        q = d.get("eccentricity", planet=n)
        t = q.time.to_value("kyr")
        e = q.data.value
        ax = axes[0,1]
        ax.plot(t, e)

    ax = axes[0,0]
    ax.legend()
    ax.set_yscale("log")
    ax.set_ylabel(r"|$a$-1au| / 1au")

    ax = axes[0,1]
    ax.set_ylabel(r"$e$")
    ax.set_xlabel(r"$t$ [kyr]")

    Npl = len(d.avail()["planets"])
    for n in range(Npl):
        q = d.get("x", planet=n)
        t = q.time.to_value("kyr")
        x = q.data.to_value("au")
        ax = axes[1,0]
        ax.plot(t, x, label=f"{n}")

        q = d.get("y", planet=n)
        t = q.time.to_value("kyr")
        y = q.data.to_value("au")
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