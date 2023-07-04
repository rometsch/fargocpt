#!/usr/bin/env python3

import numpy as np
import yaml
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import astropy.constants as const
import astropy.units as u

from create_input import get_fargo_grid, get_solution_array


def plot_field(Xi, Yi, Z, logplot=False, cmap=None, log_vfloor=1e-20, vmin=None, vmax=None):
    fig, ax = plt.subplots(dpi=150)

    if cmap is None:
        cmap = "magma"
    if logplot:
        if vmin is None:
            vmin = log_vfloor + np.min(np.abs(Z))
        if vmax is None:
            vmax = log_vfloor + np.max(np.abs(Z))
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
        Z = np.abs(Z)
    else:
        if vmin is None:
            vmin = np.min(Z)
        if vmax is None:
            vmax = np.max(Z)
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(Xi,Yi,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)

    return fig, ax, cbar

def Eint_to_Erad(Eint):
    """ Convert internal energy density to radiative energy density based on the test model. """
    mu = 2.35
    gamma = 1.4

    Rg = const.k_B / const.m_p
    cv = Rg / mu / (gamma - 1)
    cv = cv.cgs.value
    Sigma = 1 # g/cm3

    Tgas = Eint/cv/Sigma

    # assume instantaneous equilibrium
    Trad = Tgas

    sigmaR = const.sigma_sb
    aR = 4*sigmaR/const.c
    aR = aR.to_value("erg/(cm3*K4)")

    Erad = aR * Trad**4
    return Erad


if __name__ == "__main__":
    # now create an energy array with values in cgs with the energy located in one cell

    setupfile = "setup.yml"
    g = get_fargo_grid(setupfile)
    
    outdir = "output/out"


    # get final time
    with open(setupfile, "r") as infile:
        params = yaml.safe_load(infile)
        DT = float(params["DT"])
    # get initial time
    with open("test_settings.yml", "r") as infile:
        params = yaml.safe_load(infile)
        t0 = float(params["t0"])
    
    analytical_solution = get_solution_array("setup.yml", t0+DT)
    
    # get code solution from file
    code_solution = np.fromfile(f"{outdir}/snapshots/1/energy.dat", dtype=np.float64).reshape(g.Nrad, g.Naz)
    code_solution = Eint_to_Erad(code_solution)

    print("max analytical {:.4e}".format(np.max(analytical_solution)))
    print("max code {:.4e}".format(np.max(code_solution)))

    # xlim=[47,53]
    # ylim=[-3,3]
    xlim=ylim=None

    Z = analytical_solution/np.max(analytical_solution)
    Z = np.maximum(Z, 1e-10)
    fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, vmin=1e-5, vmax=1, logplot=True)
    ax.set_title("analytical")
    ax.set_xlim(xlim); ax.set_ylim(ylim)

    Z = code_solution/np.max(code_solution)
    fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, vmin=1e-5, vmax=1, logplot=True)
    ax.set_title("code")
    ax.set_xlim(xlim); ax.set_ylim(ylim)

    plt.show()