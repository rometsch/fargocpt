#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

from drift_theo import vdrift_theo

from fargocpt import Loader

def main():
    fig = plot_trajectory("../../output/tests/dust_drift/out")
    fig.savefig("trajectories.jpg", dpi=150, bbox_inches='tight')
    fig.savefig("trajectories.pdf", dpi=300, bbox_inches='tight')


def plot_trajectory(outdir):

    fig, ax = plt.subplots(figsize=(6,4))

    l = Loader(outdir)
    particles = l.particles.timeseries(["r", "stokes", "size"])


    sizes = [1e-2, 1e-1, 1, 1e1, 1e2,1e3]

    for i, p in particles.items():
        try:
            # if not i%10 == 0:
            #     continue

            t = p["time"]
            r = p["r"]
            
            
            size = p["size"][0].to("cm")
            
            if not any( [ np.isclose(size.to_value("cm"), s) for s in sizes]):
                continue
            
            stokes = p["stokes"]
            vtheo = vdrift_theo(stokes, r.to("au"))
            vtheo = vtheo.to("cm/s")
            
            # print(t, r, stokes)
            
            r = r.to("cm")
            t = t.to("s")
            rdot = (r[1:] - r[:-1])/(t[1:] - t[:-1])
            rdot = rdot.to("cm/s")
            
            Navg = -len(rdot)//10
            q = np.average(rdot[:-Navg]/vtheo[:-Navg])
            # print("stokes", stokes[-1], "\t size = ", size, "\t sim / theo =", q)
            
            label = f"s={size:.1e}, St={stokes[-1]:.1e}"
            
            ax = ax
            line, = ax.plot(t[:-1].to("yr"), -rdot, label=label)
            color = line.get_color()
            ax.plot(t.to("yr"), -vtheo, ls=":", color=color)
            
        except IndexError as e:
            print(e)

    # secondary y axis
    secax = ax.secondary_yaxis(
        'right', functions=(cmps_to_aupyr, aupyr_to_cmps))
    secax.set_ylabel(r"$-\dot{r}$ [au/yr]")


    ax.set_yscale("log")
    ax.set_ylabel(r"$-\dot{r}$ [cm/s]")
    ax.set_xlabel(r"$t$ [yr]")    
    ax.legend()


    return fig


def cmps_to_aupyr(x):
    factor = (1*u.au/u.yr).to("cm/s").value
    return x / factor


def aupyr_to_cmps(x):
    factor = (1*u.au/u.yr).to("cm/s").value
    return x * factor



if __name__ == "__main__":
    main()
