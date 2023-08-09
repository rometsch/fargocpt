#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider
import disgrid
import argparse


def main():
    opts = parser_cmdline()

    data = disgrid.Data(opts.outputdir)

    avail = data.avail()

    fig, axd = plt.subplot_mosaic([["Nbody", "Sigma"], ["slider", "slider"]], figsize=(12, 8), height_ratios=[1, 0.1])

    Nlast = avail["Nlast"]
    tlast = data.get(var="mass density", dim="2d",
                     N=Nlast).time.to_value("kyr")
    a = data.get(var="semi-major axis", planet=0)
    klast = np.argmin(np.abs(a.time.to_value("kyr") - tlast))

    # Nbody
    ax = axd["Nbody"]
    ax.set_xlabel("Time [kyr]")
    ax.set_ylabel(r"$a$ [au]")
    planet_markers = []
    planet_lines = []
    planet_fills = []
    for planet in avail["planets"]:
        a = data.get(var="semi-major axis", planet=planet)
        line, = ax.plot(a.time.to_value("kyr"), a.data.to_value(
            "au"), label=f"planet {planet}")
        planet_lines.append(line)
        e = data.get(var="eccentricity", planet=planet).data
        rmin = a.data.to_value("au") * (1 - e)
        rmax = a.data.to_value("au") * (1 + e)
        lines = ax.fill_between(a.time.to_value("kyr"), rmin, rmax,color=line.get_color(),alpha=0.4)
        planet_fills.append(lines)
        line, = ax.plot([a.time.to_value("kyr")[klast]], [a.data.to_value(
            "au")[klast]], "o", markersize=10, color=line.get_color())
        planet_markers.append(line)
    ax.legend()

    # Disk

    # Sigma
    ax = axd["Sigma"]
    field = data.get(var="mass density", dim="2d", N=avail["Nlast"])
    ax.set_xlabel("x [au]")
    ax.set_ylabel("y [au]")
    ax.set_aspect("equal")
    R, Phi = field.grid.get_meshgrid()
    X = R * np.cos(Phi)
    Y = R * np.sin(Phi)

    Z = field.data.to_value("g/cm**2")
    my_cmap = mpl.colormaps.get_cmap("viridis")
    my_norm = mpl.colors.LogNorm(vmin=np.min(Z), vmax=np.max(Z))
    pm = ax.pcolormesh(X.to_value("au"), Y.to_value("au"),
                       Z, norm=my_norm, cmap=my_cmap)
    cbar = fig.colorbar(pm, ax=ax)
    cbar.set_label(r"$\Sigma$ [g/cm$^2$]")

    fig.suptitle(f"N = {Nlast}, t = {tlast:.2e} kyr")

    #
    # add Slider
    #

    # Make a horizontal sliders to control the snapshot number.
    Nlast = int(avail["Nlast"])
    Nfirst = int(avail["Nfirst"])

    N_slider = Slider(
        ax=axd["slider"],
        label='N',
        valmin=Nfirst,
        valmax=Nlast,
        valinit=Nlast,
        valstep=1,
    )

    #
    # Add update functions for interactivity
    #

    def update_Nbody(tnow):
        for planet, marker, line, fill in zip(avail["planets"], planet_markers, planet_lines, planet_fills):
            a = data.get(var="semi-major axis", planet=planet)
            know = np.argmin(np.abs(a.time.to_value("kyr") - tnow))
            marker.set_xdata([a.time.to_value("kyr")[know]])
            marker.set_ydata([a.data.to_value("au")[know]])
            line.set_xdata(a.time.to_value("kyr"))
            line.set_ydata(a.data.to_value("au"))

            e = data.get(var="eccentricity", planet=planet).data
            rmin = a.data.to_value("au") * (1 - e)
            rmax = a.data.to_value("au") * (1 + e)
            dummy = ax.fill_between(a.time.to_value("kyr"), rmin, rmax, alpha=0)
            #create invisible dummy object to extract the vertices 
            dp = dummy.get_paths()[0]
            dummy.remove()
            #update the vertices of the PolyCollection
            fill.set_paths([dp.vertices])

    def update_Sigma(sigma):
        Z = sigma.data.to_value("g/cm**2")
        # my_nom = # you will need to scale your read data between [0, 1]
        new_data = my_norm(Z)
        new_color = my_cmap(new_data)
        pm.update({'array': new_color})

    # The function to be called anytime a slider's value changes
    def update(val, N=None):
        if N is None:
            N = int(N_slider.val)
        if N > Nlast or N < Nfirst:
            return
        sigma = data.get(var="mass density", dim="2d", N=N)
        tnow = sigma.time.to_value("kyr")
        update_Nbody(tnow)
        update_Sigma(sigma)
        fig.suptitle(f"N = {N}, t = {tnow:.2e} kyr")
        N_slider.set_val(N)
        fig.canvas.draw_idle()

    # Update the plot when mouse button is released instead of on update.
    # Otherwise the widget quickly becomes unresponsive due to exessive optool calls.
    fig.canvas.mpl_connect("button_release_event", update)

    # Update on keypress when left or right arrow key is pressed.
    def on_press(event):
        if event.key == 'right':
            update(event, N=N_slider.val+1)
        elif event.key == 'left':
            update(event, N=N_slider.val-1)
        elif event.key == 'q':
            plt.close(fig)
            exit(0)

    fig.canvas.mpl_connect('key_press_event', on_press)

    plt.show(block=False)

    # check current last N in simulation directory
    N_last_old = np.genfromtxt(f"{opts.outputdir}/snapshots/list.txt", dtype=int)[-1]

    while True:
        plt.pause(2)
        N_last_new = np.genfromtxt(f"{opts.outputdir}/snapshots/list.txt", dtype=int)[-1]
        if N_last_new > N_last_old and N_slider.val == Nlast:
            N_last_old = N_last_new
            data = disgrid.Data(opts.outputdir)
            avail = data.avail()
            Nlast = int(avail["Nlast"])
            N_slider.valmax = Nlast
            update(None, N=Nlast)


def parser_cmdline():
    parser = argparse.ArgumentParser(
        description="Plotting an overview of the simulation.")
    parser.add_argument("--outputdir", type=str,
                        default="output/out", help="Path to the data file.")
    opts = parser.parse_args()
    return opts


if __name__ == "__main__":
    main()
