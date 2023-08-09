#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider
import disgrid
import argparse


def main():
    opts = parser_cmdline()
    overview = Overview(opts.outputdir, opts.update_interval)
    overview.widget()
    overview.show(follow=True)


class Overview:

    def __init__(self, outputdir, update_interval=0.1):
        self.outputdir = outputdir
        self.update_interval = update_interval

        self.updates = []

    def show(self, follow=False):

        if not follow:
            plt.show()
        else:
            plt.show(block=False)

            # check current last N in simulation directory
            N_last_old = np.genfromtxt(
                f"{self.outputdir}/snapshots/list.txt", dtype=int)[-1]

            while True:
                # plt.pause(2)
                # time.sleep(2)
                self.fig.canvas.start_event_loop(self.update_interval)
                N_last_new = np.genfromtxt(
                    f"{self.outputdir}/snapshots/list.txt", dtype=int)[-1]
                if N_last_new > N_last_old and self.N_slider.val == self.Nlast:
                    N_last_old = N_last_new
                    self.data = disgrid.Data(self.outputdir)
                    self.Nlast = int(self.data.avail()["Nlast"])
                    self.update(None, N=self.Nlast)

    def widget(self):

        self.data = disgrid.Data(self.outputdir)
        data = self.data

        avail = data.avail()

        self.fig, self.axd = plt.subplot_mosaic(
            [["Nbody", "Sigma"], ["slider", "slider"]], figsize=(12, 8), height_ratios=[1, 0.1])

        self.Nfirst = int(avail["Nfirst"])
        self.Nlast = avail["Nlast"]
        self.Nnow = self.Nlast

        self.create_Nbody()
        self.updates.append(self.update_Nbody)

        self.create_Sigma()
        self.updates.append(self.update_Sigma)

        self.tnow = self.data.get(
            var="mass density", dim="2d", N=self.Nlast).time.to_value("kyr")
        self.update_title()

        self.create_slider()
        self.register_keys()

    def create_slider(self):
        # Make a horizontal sliders to control the snapshot number.
        self.N_slider = Slider(
            ax=self.axd["slider"],
            label='N',
            valmin=self.Nfirst,
            valmax=self.Nlast,
            valinit=self.Nlast,
            valstep=1,
        )
        # Update the plot when mouse button is released instead of on update.
        # Otherwise the widget quickly becomes unresponsive due to exessive optool calls.
        self.fig.canvas.mpl_connect("button_release_event", self.update)

    def register_keys(self):
        # Update on keypress when left or right arrow key is pressed.
        def on_press(event):
            if event.key == 'right':
                self.update(event, N=self.N_slider.val+1)
            elif event.key == 'left':
                self.update(event, N=self.N_slider.val-1)
            elif event.key == 'q':
                plt.close(self.fig)
                exit(0)
        self.fig.canvas.mpl_connect('key_press_event', on_press)

    def update_slider(self):
        sl = self.N_slider
        sl.valmax = self.Nlast
        sl.ax.set_xlim(sl.valmin, sl.valmax)

    # The function to be called anytime a slider's value changes
    def update(self, val, N=None):
        if N is None:
            N = int(self.N_slider.val)
        if N > self.Nlast or N < self.Nfirst:
            return
        sigma = self.data.get(var="mass density", dim="2d", N=N)
        self.tnow = sigma.time.to_value("kyr")
        self.Nnow = N
        for up in self.updates:
            up()
        self.fig.suptitle(f"N = {N}, t = {self.tnow:.2e} kyr")
        self.N_slider.set_val(N)
        self.update_slider()
        self.fig.canvas.draw_idle()

    def update_title(self):
        self.fig.suptitle(f"N = {self.Nnow}, t = {self.tnow:.2e} kyr")

    def update_Sigma(self):
        sigma = self.data.get(var="mass density", dim="2d", N=self.Nnow)
        Z = sigma.data.to_value("g/cm**2")
        # my_nom = # you will need to scale your read data between [0, 1]
        new_data = self.my_norm(Z)
        new_color = self.my_cmap(new_data)
        self.pm.update({'array': new_color})

    def create_Sigma(self):
        ax = self.axd["Sigma"]
        field = self.data.get(var="mass density", dim="2d", N=self.Nlast)
        ax.set_xlabel("x [au]")
        ax.set_ylabel("y [au]")
        ax.set_aspect("equal")
        R, Phi = field.grid.get_meshgrid()
        X = R * np.cos(Phi)
        Y = R * np.sin(Phi)

        Z = field.data.to_value("g/cm**2")
        self.my_cmap = mpl.colormaps.get_cmap("viridis")
        self.my_norm = mpl.colors.LogNorm(vmin=np.min(Z), vmax=np.max(Z))
        self.pm = ax.pcolormesh(X.to_value("au"), Y.to_value("au"),
                                Z, norm=self.my_norm, cmap=self.my_cmap)
        cbar = self.fig.colorbar(self.pm, ax=ax)
        cbar.set_label(r"$\Sigma$ [g/cm$^2$]")

    def update_Nbody(self):
        tnow = self.tnow
        for planet, marker, line, fill in zip(self.data.avail()["planets"], self.planet_markers, self.planet_lines, self.planet_fills):
            ax = self.axd["Nbody"]

            a = self.data.get(var="semi-major axis", planet=planet)
            know = np.argmin(np.abs(a.time.to_value("kyr") - tnow))
            marker.set_xdata([a.time.to_value("kyr")[know]])
            marker.set_ydata([a.data.to_value("au")[know]])
            line.set_xdata(a.time.to_value("kyr"))
            line.set_ydata(a.data.to_value("au"))

            e = self.data.get(var="eccentricity", planet=planet).data
            rmin = a.data.to_value("au") * (1 - e)
            rmax = a.data.to_value("au") * (1 + e)
            dummy = ax.fill_between(
                a.time.to_value("kyr"), rmin, rmax, alpha=0)
            # create invisible dummy object to extract the vertices
            dp = dummy.get_paths()[0]
            dummy.remove()
            # update the vertices of the PolyCollection
            fill.set_paths([dp.vertices])

            ax.autoscale_view()

    def create_Nbody(self):
        # Nbody
        data = self.data
        tlast = data.get(var="mass density", dim="2d",
                         N=self.Nlast).time.to_value("kyr")
        a = data.get(var="semi-major axis", planet=0)
        klast = np.argmin(np.abs(a.time.to_value("kyr") - tlast))
        avail = data.avail()

        ax = self.axd["Nbody"]
        ax.set_xlabel("Time [kyr]")
        ax.set_ylabel(r"$a$ [au]")
        self.planet_markers = []
        self.planet_lines = []
        self.planet_fills = []
        for planet in avail["planets"]:
            a = data.get(var="semi-major axis", planet=planet)
            line, = ax.plot(a.time.to_value("kyr"), a.data.to_value(
                "au"), label=f"planet {planet}")
            self.planet_lines.append(line)
            e = data.get(var="eccentricity", planet=planet).data
            rmin = a.data.to_value("au") * (1 - e)
            rmax = a.data.to_value("au") * (1 + e)
            lines = ax.fill_between(a.time.to_value(
                "kyr"), rmin, rmax, color=line.get_color(), alpha=0.4)
            self.planet_fills.append(lines)
            line, = ax.plot([a.time.to_value("kyr")[klast]], [a.data.to_value(
                "au")[klast]], "o", markersize=10, color=line.get_color())
            self.planet_markers.append(line)
        ax.legend()


def parser_cmdline():
    parser = argparse.ArgumentParser(
        description="Plotting an overview of the simulation.")
    parser.add_argument("--outputdir", type=str,
                        default="output/out", help="Path to the data file.")
    parser.add_argument("--update_interval", type=float,
                        default=0.1, help="Update interval in seconds.")
    opts = parser.parse_args()
    return opts


if __name__ == "__main__":
    main()
