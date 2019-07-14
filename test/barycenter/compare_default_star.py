#!/usr/bin/env python3
import os
from subprocess import run
import numpy as np
import argparse
import simscripts as sims
from simscripts.simcodes.fargocpt.param import Parameters
import reader_fargo
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", default=False, action="store_true", help="recalculate sims")
    parser.add_argument("-n", default=1, type=int, help="number of orbits to calculate")
    parser.add_argument("--nbody", default=False, action="store_true", help="plot nbody stuff")
    parser.add_argument("--rho", default=False, action="store_true", help="plot density colormaps")
    parser.add_argument("--outdirs", nargs=2, help="output dirs to compare")
    args = parser.parse_args()

    if (not args.nbody) and (not args.rho):
        args.nbody = True
        args.rho = True
    
    if args.r:
        params_ds = Parameters("gamma-ceph-defaultstar.par")
        params_ds.set("Ntot", args.n*10)
        
        params_nods = Parameters("gamma-ceph-nodefaultstar.par")
        params_nods.set("Ntot", args.n*10)

        run_sims()

    if args.outdirs:
        outdir1 = args.outdirs[0]
        outdir2 = args.outdirs[1]
    else:
        outdir1 = "out/comp-defaultstar/"
        outdir2 = "out/comp-nodefaultstar/"
    # default star data
    reader_ds = reader_fargo.Reader(outdir1)
    # no default star data
    reader_nods = reader_fargo.Reader(outdir2)

    if args.n:
        Nrho = args.n
    else:
        Nrho = len(reader_ds.outputTimes)-1
    print("Comparing output number {}".format(Nrho))
    compare_binary_output(os.path.join(outdir1,"gasdens{}.dat".format(Nrho)),
                          os.path.join(outdir2,"gasdens{}.dat".format(Nrho)))
        
    # default star data
    reader_ds = reader_fargo.Reader(outdir1)
    # no default star data
    reader_nods = reader_fargo.Reader(outdir2)

    if args.nbody:
        plot_nbody(reader_ds, reader_nods)
    if args.rho:
        plot_rho(reader_ds, reader_nods, n=Nrho)

def plot_nbody(reader1, reader2):

    params1 = get_param(reader1)
    params2 = get_param(reader2)
    
    fig, axes = plt.subplots(2,2,constrained_layout=True); axes = axes.ravel()
    varnames = ["semi-major axis", "eccentricity", "x", "y"]
    npl = 1
    
    for k, vn in enumerate(varnames):

        ax = axes[k]
        varname = "semi-major axis"
        time, val = reader1.getScalar(vn, get_nsecondary(params1))
        ax.plot(time/4.956246157/100, val, label="ds")
        
        time, val = reader2.getScalar(vn, get_nsecondary(params2))
        ax.plot(time/4.956246157/100, val, label="nods")
        
        ax.set_ylabel(vn)
        ax.set_xlabel("time in binary orbits")

        ax.grid(color="grey", alpha=0.3)
        if k == 0:
            ax.legend()

def get_param(reader):
    for fname in os.listdir(reader.dataDir):
        if os.path.splitext(fname)[1] == ".par":
            return Parameters(os.path.join(reader.dataDir, fname))

def get_nsecondary(param):
    try:
        if param.get("NoDefaultStar") == "YES":
            return 2
        else:
            return 1
    except KeyError:
        return 1
            
def plot_rho(reader1, reader2, n):

    params1 = get_param(reader1)
    params2 = get_param(reader2)
    
    Nrho = n
    
    fig, axes = plt.subplots(1,3, constrained_layout=True)
    rho_ds = reader1.load2d(Nrho, "gasdens{}.dat", u.g/u.cm**2)
    rho_nods = reader2.load2d(Nrho, "gasdens{}.dat", u.g/u.cm**2)
    X = reader1.X.to("au")
    Y = reader1.Y.to("au")

    delta_rho = np.abs(rho_ds - rho_nods)/np.max(np.abs(rho_ds))

    vmax = np.max( [np.max(rho_ds.value), np.max(rho_nods.value)])
    vmin = np.min( [np.min(rho_ds.value), np.min(rho_nods.value)])

    cmap = "coolwarm"

    norm=colors.LogNorm(vmin=vmin, vmax=vmax)
    
    ax = axes[0]
    im = ax.pcolormesh(X,Y, rho_ds.value, vmax = vmax, cmap = cmap, norm=norm)
    ax.set_title("default star")
    add_colorbar(fig, ax, im, 'rho [g/cm2]')
    # add line pointing to planet
    nsecondary = get_nsecondary(params1)
    t,x = reader1.getScalar("x", n=nsecondary, frame=Nrho)
    t,y = reader1.getScalar("y", n=nsecondary, frame=Nrho)
    x = x.to("au").value
    y = y.to("au").value
    r = np.sqrt(x**2 + y**2)/np.max(X.value)
    ax.plot([0.0, x/r], [0.0, y/r], linewidth=0.3, color="green")

    ax = axes[1]
    im = ax.pcolormesh(X,Y, rho_nods.value, vmax = vmax, cmap = cmap, norm=norm)
    ax.set_title("no default star")
    add_colorbar(fig, ax, im, "rho [g/cm2]")
    # add line pointing to planet
    nsecondary = get_nsecondary(params2)
    t,x = reader2.getScalar("x", n=nsecondary, frame=Nrho)
    t,y = reader2.getScalar("y", n=nsecondary, frame=Nrho)
    x = x.to("au").value
    y = y.to("au").value
    r = np.sqrt(x**2 + y**2)/np.max(X.value)
    ax.plot([0.0, x/r], [0.0, y/r], linewidth=0.3, color="green")

    
    ax = axes[2]
    im = ax.pcolormesh(X,Y,np.log(delta_rho), vmin=-10, vmax=0, cmap = cmap)
    ax.set_title("log( delta rho rel)")
    add_colorbar(fig, ax, im, "log(delta rho/max rho)")
    
    for ax in axes:
        ax.set_aspect("equal")

    fig.suptitle("surface density after {} orbits".format(Nrho/10))
        
    plt.show()


def compare_binary_output(file1, file2):
    run(["../../Tools/compare_binary_output.py", file1, file2])
    
def run_sims():
    for configfile in ["gamma-ceph-defaultstar.par", "gamma-ceph-nodefaultstar.par"]:
        run(["mpirun", "-n", "4", "../../fargo", configfile])

        
def add_colorbar(fig, ax, im, label):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', size='5%', pad=0.05)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    ax.xaxis.set_major_locator(plt.NullLocator())
    cbar.set_label(label)
        
if __name__=="__main__":
    main()
