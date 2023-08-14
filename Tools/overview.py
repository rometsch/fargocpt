#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider
try:
    import disgrid
except ImportError:
    print("Dependency 'disgrid' not found. Please install it via 'pip install git+https://github.com/rometsch/disgrid'")
import argparse


def main():
    opts = parser_cmdline()
    vars=["0:Nbody", "2:mass density", "2:velocity azimuthal", "0:eccentricity", "0:mass"]
    overview = Overview(opts.outputdir, opts.follow, vars=vars)
    overview.widget()
    overview.show()


class Plot2D:

    def __init__(self, ax, simd, var):
        self.ax = ax
        self.fig = ax.figure
        self.simd = simd
        self.var = var
        self.type = "2"

    def update(self, Nnow, tnow):
        sigma = self.simd.get(var=self.var, dim="2d", N=Nnow)
        Z = sigma.data.cgs.value
        # my_nom = # you will need to scale your read data between [0, 1]
        new_data = self.my_norm(Z)
        new_color = self.my_cmap(new_data)
        self.pm.update({'array': new_color})

    def create(self):
        ax = self.ax
        Nlast = self.simd.avail()["Nlast"]
        field = self.simd.get(var=self.var, dim="2d", N=Nlast)
        ax.set_xlabel("x [au]")
        ax.set_ylabel("y [au]")
        ax.set_aspect("equal")
        R, Phi = field.grid.get_meshgrid()
        X = R * np.cos(Phi)
        Y = R * np.sin(Phi)

        Z = field.data.cgs.value
        vmin=np.min(Z)
        vmax=np.max(Z)
        if vmin >= 0:
            self.my_cmap = mpl.colormaps.get_cmap("viridis")
            self.my_norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            vabs = max(np.abs(vmin), vmax)
            vmin = -vabs
            vmax = vabs
            self.my_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            self.my_cmap = mpl.colormaps.get_cmap("bwr")
        self.pm = ax.pcolormesh(X.to_value("au"), Y.to_value("au"),
                                Z, norm=self.my_norm, cmap=self.my_cmap)
        if vmin != vmax:
            self.cbar = self.fig.colorbar(self.pm, ax=ax)
            self.cbar.set_label(f"{self.var} [{field.data.cgs.unit}]")

class PlotNbody:

    def __init__(self, ax, simd, var):
        self.ax = ax
        self.simd = simd
        self.var = var
        self.type = "0"

    def update(self, Nnow, tnow):
        for planet, line, fill in zip(self.simd.avail()["planets"], self.planet_lines, self.planet_fills):
            ax = self.ax
            a = self.simd.get(var="semi-major axis", planet=planet)
            line.set_xdata(a.time.to_value("kyr"))
            line.set_ydata(a.data.to_value("au"))

            e = self.simd.get(var="eccentricity", planet=planet).data
            rmin = a.data.to_value("au") * (1 - e)
            rmax = a.data.to_value("au") * (1 + e)
            dummy = ax.fill_between(
                a.time.to_value("kyr"), rmin, rmax, alpha=0)
            # create invisible dummy object to extract the vertices
            dp = dummy.get_paths()[0]
            dummy.remove()
            # update the vertices of the PolyCollection
            fill.set_paths([dp.vertices])


        self.timemarker.set_xdata([tnow, tnow])
        ax.autoscale_view()

    def create(self):
        # Nbody
        data = self.simd
        tlast = data.loader.snapshot_times[-1].to_value("kyr")
        a = data.get(var="semi-major axis", planet=0)
        avail = data.avail()

        ax = self.ax
        ax.set_xlabel("Time [kyr]")
        ax.set_ylabel(r"$a$ [au]")
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
        ax.legend()

        self.timemarker = ax.axvline(x=tlast, color="k", linestyle="-", alpha=0.5)



class PlotScalar:

    def __init__(self, ax, simd, vars):
        self.ax = ax
        self.simd = simd
        self.vars = vars
        self.type = "0"

    def update(self, Nnow, tnow):
        ax = self.ax
        for var, line in zip(self.vars, self.lines):
            x = self.simd.get(var=var, dim="scalar")

            line.set_xdata(x.time.to_value("kyr"))
            line.set_ydata(x.data.cgs.value)

        self.timemarker.set_xdata([tnow, tnow])
        ax.autoscale_view()



    def create(self):
        # Nbody
        data = self.simd
        tlast = data.loader.snapshot_times[-1].to_value("kyr")
        x = data.get(var="mass", dim="scalar")
        klast = np.argmin(np.abs(x.time.to_value("kyr") - tlast))

        ax = self.ax
        ax.set_xlabel("Time [kyr]")
        ax.set_ylabel(r"value [cgs]")

        self.lines = []

        for var in self.vars:
            ax = self.ax
            x = self.simd.get(var=var, dim="scalar")
            line, = ax.plot(x.time.to_value("kyr"), x.data.cgs, 
                            label=f"{var}")
            self.lines.append(line)

        ax.legend()

        self.timemarker = ax.axvline(x=tlast, color="k", linestyle="-", alpha=0.5)


class Overview:

    def __init__(self, outputdir, update_interval=0.0, vars=["0:Nbody", "2:mass density", "2:velocity azimuthal"]):
        self.outputdir = outputdir
        self.update_interval = update_interval
        self.vars = [k.split(":")[-1] for k in vars]
        self.plot_types = [k.split(":")[0] for k in vars]

    def show(self, follow=None):

        if follow is None:
            follow = self.update_interval

        if follow == 0.0:
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

        Nplots = len(self.vars)
        Ncols = int(np.ceil(np.sqrt(Nplots)))
        Nrows = int(np.ceil(Nplots / Ncols))

        mosaic_def =  []
        for k in range(Nrows):
            rowdef = self.vars[k*Ncols:(k+1)*Ncols]
            if len(rowdef) < Ncols:
                rowdef += ["."] * (Ncols - len(rowdef))
            mosaic_def.append(rowdef)
        
        height_ratios = [1] * Nrows + [0.1]

        mosaic_def.append(["slider"] * Ncols)

        self.fig, self.axd = plt.subplot_mosaic(
            mosaic_def, figsize=(12, 8), height_ratios=height_ratios)

        self.Nfirst = int(avail["Nfirst"])
        self.Nlast = avail["Nlast"]
        self.Nnow = self.Nlast


        self.plots = {}

        for v, t in zip(self.vars, self.plot_types):
            if v == "Nbody":
                self.plots[v] = PlotNbody(self.axd[v], data, v)
            elif t == "2":
                self.plots[v] = Plot2D(self.axd[v], data, v)
            elif t == "0":
                self.plots[v] = PlotScalar(self.axd[v], data, [s.strip() for s in v.split(',')])
            else:
                NotImplementedError(f"Plot type {t} not implemented for variable {v}")

        for plot in self.plots.values():
            plot.create()

        self.tnow = self.data.get(
            var="mass density", dim="2d", N=self.Nlast).time.to_value("kyr")
        self.update_title(self.Nnow, self.tnow)

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
        for name, plot in self.plots.items():
            plot.simd = self.data
            plot.update(self.Nnow, self.tnow)
        self.fig.suptitle(f"N = {N}, t = {self.tnow:.2e} kyr")
        self.N_slider.set_val(N)
        self.update_slider()
        self.fig.canvas.draw_idle()

    def update_title(self, Nnow, tnow):
        self.fig.suptitle(f"N = {Nnow}, t = {tnow:.2e} kyr")


def parser_cmdline():
    parser = argparse.ArgumentParser(
        description="Plotting an overview of the simulation.")
    parser.add_argument("outputdir", type=str,
                        default="output/out", help="Path to the data file.")
    parser.add_argument("-f", "--follow", type=float,
                        default=0, help="Update interval in seconds. Disable if 0.")
    opts = parser.parse_args()
    return opts


if __name__ == "__main__":
    main()
