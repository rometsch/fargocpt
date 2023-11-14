#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider
from fargocpt import Loader
import argparse


default_vars = ["0:Nbody", "2:Sigma:rphi", "2:vazi:rphi", "2:vrad:rphi", "0:mass"]

def main():
    opts = parser_cmdline()
    overview = Overview(opts.outputdir, opts.follow, vars=opts.vars, start=opts.start)
    overview.create()
    overview.show()
    plt.show()


class Plot2D:

    def __init__(self, ax, loader, var, spec, start=None):
        self.ax = ax
        self.fig = ax.figure
        self.loader = loader
        self.var = var
        self.type = "2"
        self.spec = spec
        self.start = start
        parts = spec.split(":")
        self.rphi = len(parts) > 2 and "rphi" in parts[2]
        self.xy = not self.rphi
        self.plot_relative = len(parts) > 2 and "rel" in parts[2]
        self.plot_diff = len(parts) > 2 and "diff" in parts[2]
        self.label_modifier = ""
        if self.plot_relative:
            self.label_modifier += "rel"
        if self.plot_diff:
            self.label_modifier += "diff"

    def update(self, Nnow, tnow):
        l = self.loader
        Z = l.gas.vars2D.get(self.var, Nnow, grid=False).cgs.value
        if self.plot_relative:
            Z = Z/self.Z0 -1
        if self.plot_diff:
            Z = Z - self.Z0
        if Z.shape[0] == self.X.shape[0]:
            Z = Z[:-1,:]
        # my_nom = # you will need to scale your read data between [0, 1]
        new_data = self.my_norm(Z)
        new_color = self.my_cmap(new_data)
        self.pm.update({'array': new_color})

    def create(self):
        l = self.loader
        ax = self.ax
        if self.start is None:
            N = l.snapshots[-1]
        else:
            N = self.start

        if self.rphi:
            ax.set_xlabel("r [au]")
            ax.set_ylabel("phi [rad]")
        else:
            ax.set_xlabel("x [au]")
            ax.set_ylabel("y [au]")
            ax.set_aspect("equal")

        R, Phi, Z = l.gas.vars2D.get(self.var, N, grid_for_plot=True)
        dataunit = Z.cgs.unit
        Z = Z.cgs.value
        R = R.to_value("au")

        if self.rphi:
            X = R
            Y = Phi
        else:
            X = R * np.cos(Phi)
            Y = R * np.sin(Phi)

        self.X = X
        self.Y = Y

        if self.plot_relative or self.plot_diff:
            try:
                vals = l.gas.vars2D.get(self.var, "reference", grid=False)
            except KeyError:
                vals = l.gas.vars2D.get(self.var, 0, grid=False)
            self.Z0 = vals.cgs.value
        if self.plot_relative:
            Z = Z/self.Z0 -1
        if self.plot_diff:
            Z = Z - self.Z0

        if Z.shape[0] == self.X.shape[0]:
            Z = Z[:-1,:]

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
        self.pm = ax.pcolormesh(X, Y,
                                Z, norm=self.my_norm, cmap=self.my_cmap)
        if vmin != vmax:
            self.cbar = self.fig.colorbar(self.pm, ax=ax)
            self.cbar.set_label(f"{self.var} {self.label_modifier} [{dataunit}]")


class Plot1D:

    def __init__(self, ax, loader, vars, spec, start=None):
        self.ax = ax
        self.fig = ax.figure
        self.loader = loader
        self.vars = vars
        self.type = "1"
        self.lines = []
        self.spec = spec
        self.start = start
        parts = spec.split(":")
        self.plot_diff = len(parts) > 2 and "diff" in parts[2]
        self.plot_rel = len(parts) > 2 and "rel" in parts[2]
        if self.plot_diff and self.plot_rel:
            raise ValueError("spec must not contain both 'diff' and 'rel' as modifier, Use one, e.g. '1:mass density:diff'")
        self.plot_minmax = len(parts) > 2 and "minmax" in parts[2]
        self.fills = dict()
        self.label_modifier = ""
        if self.plot_rel:
            self.label_modifier += "rel"
        if self.plot_diff:
            self.label_modifier += "diff"

    def update(self, Nnow, tnow):
        ax = self.ax

        for var, line in zip(self.vars, self.lines):
            y = self.loader.gas.vars2D.avg(var, Nnow, grid=False).cgs.value
            if self.plot_rel:
                y = y/self.y0s[var] - 1
            elif self.plot_diff:
                y = y - self.y0s[var]
            line.set_ydata(y)

            if self.plot_minmax:
                ymin = self.loader.gas.vars2D.min(var, Nnow, grid=False).cgs.value
                r, ymax = self.loader.gas.vars2D.max(var, Nnow, grid=True).cgs.value
                ymax = ymax.cgs.value
                x = r.to_value("au")
                if self.plot_rel:
                    ymin = ymin/self.y0s[var] - 1
                    ymax = ymax/self.y0s[var] - 1
                dummy = ax.fill_between(x, ymin, ymax, alpha=0)
                # create invisible dummy object to extract the vertices
                dp = dummy.get_paths()[0]
                dummy.remove()
                # update the vertices of the PolyCollection
                self.fills[var].set_paths([dp.vertices])

        self.ax.autoscale_view()

    def create(self):
        if self.start is None:
            Nnow = self.loader.snapshots[-1]
        else:
            Nnow = self.start

        ax = self.ax
        ax.set_xlabel("r [au]")

        if self.plot_rel:
            label = "y/y0 - 1"
        elif self.plot_diff:
            label = "y - y0"
        else:
            label = "value [cgs]"
        ax.set_ylabel(label)


        if self.plot_rel or self.plot_diff:
            self.y0s = {}
            for var in self.vars:
                try:
                    y0 = self.loader.gas.vars2D.avg(var, "reference", grid=False)
                except KeyError:
                    y0 = self.loader.gas.vars2D.avg(var, 0, grid=False)
                self.y0s[var] = y0.cgs.value

            if self.plot_minmax:
                self.ymax0s = {}
                self.ymin0s = {}
                for var in self.vars:
                    try:
                        ymax0 = self.loader.gas.vars2D.max(var, "reference", grid=False)
                        ymin0 = self.loader.gas.vars2D.min(var, "reference", grid=False)
                    except KeyError:
                        ymax0 = self.loader.gas.vars2D.max(var, 0, grid=False)
                        ymin0 = self.loader.gas.vars2D.min(var, 0, grid=False)
                    self.ymax0s[var] = ymax0.cgs.value
                    self.ymin0s[var] = ymin0.cgs.value


        for var in self.vars:
            ax = self.ax
            x, y = self.loader.gas.vars2D.avg(var, Nnow)
            y = y.cgs.value
            x = x.to_value("au")
            if self.plot_rel:
                y = y/self.y0s[var] - 1
            elif self.plot_diff:
                y = y - self.y0s[var]
            line, = ax.plot(x, y, label=f"{var}")
            self.lines.append(line)
            if self.plot_minmax:
                ymin = self.loader.gas.vars2D.min(var, Nnow, grid=False).cgs.value
                ymax = self.loader.gas.vars2D.max(var, Nnow, grid=False).cgs.value
                if self.plot_rel:
                    ymin = ymin/self.y0s[var] - 1
                    ymax = ymax/self.y0s[var] - 1
                self.fills[var] = ax.fill_between(x, ymin, ymax, alpha=0.4, color=line.get_color())

        ax.legend()


class PlotNbody:

    def __init__(self, ax, loader, var, spec, start=None):
        self.ax = ax
        self.loader = loader
        self.var = var
        self.type = "0"
        self.spec = spec
        self.start = start

    def update(self, Nnow, tnow):
        for nb, line, fill in zip(self.loader.nbody, self.planet_lines, self.planet_fills):
            ax = self.ax
            t = nb.time.to_value("kyr")
            a = nb.semi_major_axis.to_value("au")
            line.set_xdata(t)
            line.set_ydata(a)

            e = nb.eccentricity.value
            rmin = a * (1 - e)
            rmax = a * (1 + e)
            dummy = ax.fill_between(t, rmin, rmax, alpha=0)
            # create invisible dummy object to extract the vertices
            dp = dummy.get_paths()[0]
            dummy.remove()
            # update the vertices of the PolyCollection
            fill.set_paths([dp.vertices])


        self.timemarker.set_xdata([tnow, tnow])
        ax.autoscale_view()

    def create(self):
        # Nbody
        l = self.loader
        if self.start is None:
            N = l.snapshots[-1]
        else:
            N = self.start
        tlast = l.snapshot_time[N].to_value("kyr")

        ax = self.ax
        ax.set_xlabel("Time [kyr]")
        ax.set_ylabel(r"$a$ [au]")
        self.planet_lines = []
        self.planet_fills = []
        for nb in l.nbody:
            t = nb.time.to_value("kyr")
            a = nb.semi_major_axis.to_value("au")

            line, = ax.plot(t, a, label=f"planet {nb.id}")
            self.planet_lines.append(line)
            e = nb.eccentricity.value
            rmin = a * (1 - e)
            rmax = a * (1 + e)
            lines = ax.fill_between(t, rmin, rmax, color=line.get_color(), alpha=0.4)
            self.planet_fills.append(lines)
        ax.legend()

        self.timemarker = ax.axvline(x=tlast, color="k", linestyle="-", alpha=0.5)



class PlotScalar:

    def __init__(self, ax, loader: Loader, vars, spec, start=None):
        self.ax = ax
        self.loader = loader
        self.vars = vars
        self.type = "0"
        self.spec = spec
        self.start = start

    def update(self, Nnow, tnow):
        ax = self.ax
        for var, line in zip(self.vars, self.lines):
            x = self.loader.get(var=var, dim="scalar")
            x = self.loader.gas.scalars.get(var).cgs.value
            t = self.loader.gas.scalars.time.to_value("kyr")

            line.set_xdata(t)
            line.set_ydata(x)

        self.timemarker.set_xdata([tnow, tnow])
        ax.autoscale_view()



    def create(self):
        # Nbody
        l = self.loader
        if self.start is None:
            N = l.snapshots[-1]
        else:
            N = self.start
        tlast = l.snapshot_time(N).to_value("kyr")

        ax = self.ax
        ax.set_xlabel("Time [kyr]")
        ax.set_ylabel(r"value [cgs]")

        self.lines = []

        for var in self.vars:
            ax = self.ax
            x = l.gas.scalars.get(var).cgs.value
            t = x.time.to_value("kyr")
            line, = ax.plot(t, x, 
                            label=f"{var}")
            self.lines.append(line)

        ax.legend()

        self.timemarker = ax.axvline(x=tlast, color="k", linestyle="-", alpha=0.5)


class Overview:

    def __init__(self, outputdir, update_interval=0.0, vars=["0:Nbody", "2:Sigma:rphi", "2:vazi"], start=None, figsize=(8,6), dpi=150):
        self.outputdir = outputdir
        self.update_interval = update_interval
        self.vars = [k.split(":")[1] for k in vars]
        self.keys = vars
        self.plot_types = [k.split(":")[0] for k in vars]
        self._is_widget_created = False
        self.start = start
        self.figsize = figsize
        self.dpi = dpi

    def show(self, follow=None):

        if not self._is_widget_created:
            self.create(self.dpi)

        if follow is None or follow:
            follow = self.update_interval

        if follow == 0.0:
            self.fig.show()
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
                    self.loader = Loader(self.outputdir)
                    self.Nlast = int(self.loader.snapshots[-1])
                    self.update(None, N=self.Nlast)

    def create(self, dpi=150):

        self.loader = Loader(self.outputdir)
        ldr = self.loader

        Nplots = len(self.vars)
        Ncols = int(np.ceil(np.sqrt(Nplots)))
        Nrows = int(np.ceil(Nplots / Ncols))

        mosaic_def =  []
        for k in range(Nrows):
            rowdef = self.keys[k*Ncols:(k+1)*Ncols]
            if len(rowdef) < Ncols:
                rowdef += ["."] * (Ncols - len(rowdef))
            mosaic_def.append(rowdef)
        
        height_ratios = [1] * Nrows + [0.1]

        mosaic_def.append(["slider"] * Ncols)

        self.fig, self.axd = plt.subplot_mosaic(
            mosaic_def, figsize=self.figsize, height_ratios=height_ratios, dpi=dpi)

        self.Nfirst = ldr.snapshots[0]
        self.Nlast = ldr.snapshots[-1]
        if self.start is None:
            self.Nnow = self.Nlast
        else:
            self.Nnow = self.start


        self.plots = {}

        for v, t, k in zip(self.vars, self.plot_types, self.keys):
            if v == "Nbody":
                self.plots[k] = PlotNbody(self.axd[k], ldr, v, k, start=self.start)
            elif t == "2":
                self.plots[k] = Plot2D(self.axd[k], ldr, v, k, start=self.start)
            elif t == "0":
                self.plots[k] = PlotScalar(self.axd[k], ldr, [s.strip() for s in v.split(',')], k, start=self.start)
            elif t == "1":
                self.plots[k] = Plot1D(self.axd[k], ldr, [s.strip() for s in v.split(',')], k, start=self.start)
            else:
                NotImplementedError(f"Plot type {t} not implemented for variable {v}")

        for plot in self.plots.values():
            plot.create()

        self.tnow = ldr.snapshot_time[self.Nnow].to_value("kyr")
        self.update_title(self.Nnow, self.tnow)

        self.create_slider()
        self.register_keys()


        self.fig.tight_layout()

        self._is_widget_created = True

    def create_slider(self):
        if self.start is None:
            Ninit = self.Nlast
        else:
            Ninit = self.start

        # Make a horizontal sliders to control the snapshot number.
        self.N_slider = Slider(
            ax=self.axd["slider"],
            label='N',
            valmin=self.Nfirst,
            valmax=max(self.Nlast, self.Nfirst+1),
            valinit=Ninit,
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
        try:
            self.tnow = self.loader.snapshot_time[N].to_value("kyr")
            self.Nnow = N
            for name, plot in self.plots.items():
                plot.loader = self.loader
                plot.update(self.Nnow, self.tnow)
            self.fig.suptitle(f"N = {N}, t = {self.tnow:.2e} kyr")
            self.N_slider.set_val(N)
            self.update_slider()
            self.fig.canvas.draw_idle()
        except (TypeError, FileNotFoundError, IndexError, KeyError):
            pass

    def update_title(self, Nnow, tnow):
        self.fig.suptitle(f"N = {Nnow}, t = {tnow:.2e} kyr")


def parser_cmdline():
    parser = argparse.ArgumentParser(
        description="Plotting an overview of the simulation.")
    parser.add_argument("outputdir", type=str,
                        default="output/out", help="Path to the data file.")
    parser.add_argument("-f", "--follow", type=float,
                        default=0, help="Update interval in seconds. Disable if 0.")
    parser.add_argument("--vars", type=str, default=default_vars, nargs="+", help="Additional variables to be plotted. Specify them like '<dim=0,1,2>:<variable name>")
    parser.add_argument("-s", "--start", type=int, help="Start at this snapshot id.")
    opts = parser.parse_args()
    return opts


if __name__ == "__main__":
    main()
