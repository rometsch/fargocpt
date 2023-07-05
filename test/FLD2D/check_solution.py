#!/usr/bin/env python3

import numpy as np
import yaml
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import astropy.constants as const
import astropy.units as u
import os

from argparse import ArgumentParser

from create_input import get_fargo_grid, get_solution_array


def plot_field(Xi, Yi, Z, scaling="linear", cmap=None, log_vfloor=1e-20, vmin=None, vmax=None, linthresh=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(dpi=150)
    else:
        fig = ax.get_figure()

    if cmap is None:
        cmap = "magma"
    if scaling == "log":
        if vmin is None:
            vmin = log_vfloor + np.min(np.abs(Z))
        if vmax is None:
            vmax = log_vfloor + np.max(np.abs(Z))
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
        Z = np.abs(Z)
    elif scaling == "symlog":
        if vmax is None:
            vmax = np.max(np.abs(Z))
        if linthresh is None:
            if vmin is None:
                linthresh = 1e-2*vmax
            else:
                linthresh = vmin
        else:
            linthresh = linthresh
        norm = mplcolors.SymLogNorm(vmin=-vmax, vmax=vmax, linthresh=linthresh)
        Z = np.abs(Z)
    elif scaling == "linear":
        if vmin is None:
            vmin = np.min(Z)
        if vmax is None:
            vmax = np.max(Z)
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError(f"Scaling '{scaling}' not known.")

    pcm = ax.pcolormesh(Xi,Yi,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax, location="right")

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

def setup_hash():
    import hashlib
    with open("test_settings.yml", "r") as infile:
        test_settings = infile.read()
    with open("setup.yml", "r") as infile:
        params = yaml.safe_load(infile)
        DT = str(params["DT"])
        Nrad = str(params["Nrad"])
        Naz = str(params["Nsec"])
        vname = str(params["RadiativeDiffusionVariable"])

    s = test_settings + DT + Nrad + Naz + vname
    h = hashlib.sha256(s.encode("utf-8"))
    return h.hexdigest()
    

if __name__ == "__main__":
    # now create an energy array with values in cgs with the energy located in one cell

    parser = ArgumentParser()
    parser.add_argument("--solve", action="store_true", help="Solve the linear system with an exact method.")
    parser.add_argument("--exactmat", action="store_true", help="Also display exact matrix inversion solution.")
    parser.add_argument("--diff", action="store_true", help="Show differences.")
    opts = parser.parse_args()

    if opts.solve:
        opts.exactmat = True

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
    vmax = 1e5
    vmin = 1 

    ncols = 3 if opts.solve else 2
    nrows = 2 if opts.diff else 1
    figwidth = 12 if opts.solve else 8
    figheight = 6.5 if opts.diff else 5
    fig, axsfull = plt.subplots(ncols=ncols, nrows=nrows, dpi=150, sharex="all", sharey="all", figsize=(figwidth,figheight))
    
    if opts.diff:
        axsv = axsfull[0,:]
    else:
        axsv = axsfull


    Z = analytical_solution
    # Z = Z/np.max(Z)
    fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, scaling="log", ax=axsv[0], vmin=vmin, vmax=vmax)
    ax.set_title("analytical")
    cbar.set_label(r"$E_\mathrm{rad}$ [cgs]")
    ax.set_xlim(xlim); ax.set_ylim(ylim)

    Z = code_solution
    # Z = Z/np.max(Z)
    fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, scaling="log", ax=axsv[1], vmin=vmin, vmax=vmax)
    ax.set_title("code")
    cbar.set_label(r"$E_\mathrm{rad}$ [cgs]")
    ax.set_xlim(xlim); ax.set_ylim(ylim)

    if opts.exactmat or opts.solve:
        # check if there is a solution already saved
        # get a hash for the initial condition and DT
        h = setup_hash()
        filename = f"matexact_{h}.npy"
        if not os.path.exists(filename):
            if not opts.solve:
                print("Exact matrix solution not yet computed. Please run './check_solution -- solve' first.")
            else:
                print("Computing exact solution of linear system.")
                from solve import solve_linear_system
                initial_condition = get_solution_array("setup.yml", t0)
                Nr, Naz = initial_condition.shape
                matrix_inversion_solution = solve_linear_system(outdir, initial_condition.flatten(), Nr, Naz).reshape(Nr, Naz)
                np.save(filename, matrix_inversion_solution)
        if os.path.exists(filename):
            matrix_inversion_solution = np.load(filename)
            print("max mat exact {:.4e}".format(np.max(matrix_inversion_solution)))
            Z = matrix_inversion_solution
            # Z = Z/np.max(Z)
            fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, scaling="log", ax=axsv[2], vmin=vmin, vmax=vmax)
            ax.set_title("mat exact")
            cbar.set_label(r"$E_\mathrm{rad}$ [cgs]")
            ax.set_xlim(xlim); ax.set_ylim(ylim)

    #
    # Show differences
    #

    if opts.diff:
        axsd = axsfull[1,:]
        
        Z = code_solution - analytical_solution
        diffmax = np.max(np.abs(Z))
        fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, scaling="linear", ax=axsd[0], vmin=-diffmax, vmax=diffmax, cmap="bwr")
        ax.set_title("diff code - analytical")
        cbar.set_label(r"$E_\mathrm{rad}$ [cgs]")
        ax.set_xlim(xlim); ax.set_ylim(ylim)


    if opts.diff and opts.solve:

        Z = code_solution - matrix_inversion_solution
        diffmax = np.max(np.abs(Z))
        fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, scaling="linear", ax=axsd[1], vmin=-diffmax, vmax=diffmax, cmap="bwr")
        ax.set_title("diff code - exact mat")
        cbar.set_label(r"$E_\mathrm{rad}$ [cgs]")
        ax.set_xlim(xlim); ax.set_ylim(ylim)

        Z = matrix_inversion_solution - analytical_solution
        diffmax = np.max(np.abs(Z))
        fig, ax, cbar = plot_field(g.Xi, g.Yi, Z, scaling="linear", ax=axsd[2], vmin=-diffmax, vmax=diffmax, cmap="bwr")
        ax.set_title("exact mat - analytical")
        cbar.set_label(r"$E_\mathrm{rad}$ [cgs]")
        ax.set_xlim(xlim); ax.set_ylim(ylim)



    fig.suptitle("Comparison between analytical solution and matrix solver")

    plt.show()