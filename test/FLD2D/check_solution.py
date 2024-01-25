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

basedir = "../../output/tests/FLD2D"

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

    pcm = ax.pcolormesh(Xi, Yi, Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax, location="right")

    return fig, ax, cbar


def setup_hash():
    import hashlib
    with open("test_settings.yml", "r") as infile:
        test_settings = infile.read()
    with open("setup.yml", "r") as infile:
        params = yaml.safe_load(infile)
        Nrad = str(params["Nrad"])
        Naz = str(params["Nsec"])

    s = test_settings + Nrad + Naz
    h = hashlib.sha256(s.encode("utf-8"))
    return h.hexdigest()


if __name__ == "__main__":
    # now create an energy array with values in cgs with the energy located in one cell

    parser = ArgumentParser()
    parser.add_argument("--test", action="store_true", help="Run a pass fail test which calculates the integrated absolute difference and compares to a threshold in test_settings.yml. This prints a success or fail message.")
    parser.add_argument("--solve", action="store_true",
                        help="Solve the linear system with an exact method.")
    parser.add_argument("--exactmat", action="store_true",
                        help="Also display exact matrix inversion solution.")
    parser.add_argument("--diff", action="store_true",
                        help="Show differences.")
    parser.add_argument("--outfile", type=str, help="Filename for the plot output.")
    parser.add_argument("--verbose", action="store_true", help="Print more information.")
    opts = parser.parse_args()

    if opts.solve:
        opts.exactmat = True

    setupfile = "setup.yml"
    g = get_fargo_grid(setupfile)

    outdir = basedir + "/out"

    # get final time
    with open(setupfile, "r") as infile:
        params = yaml.safe_load(infile)
    # get initial time
    with open("test_settings.yml", "r") as infile:
        params = yaml.safe_load(infile)
    offset = float(params["offset"])

    tfinal = float(params["tfinal"])

    analytical_solution = get_solution_array("setup.yml", tfinal)

    # get code solution from file
    code_solution = np.fromfile(
        f"{outdir}/f_FLD2Dtest_output.dat", dtype=np.float64).reshape(g.Nrad, g.Naz)

    # calculate integral of functions and of differences by summing up the function value times the area of each cell
    integral_analytical = np.sum((analytical_solution - offset)*g.A)
    integral_code = np.sum((code_solution-offset)*g.A)
    integral_diff = np.sum((code_solution-analytical_solution)*g.A)
    integral_absdiff = np.sum(np.abs(code_solution-analytical_solution)*g.A)

    if opts.test:
        threshold = float(params["threshold"])
        test_name = "FLD2D"
        with open("test.log", "w") as f:
            from datetime import datetime
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{current_time}", file=f)
            print(f"integral diff = {integral_diff}, theshold = {threshold}", file=f)

        if integral_absdiff < threshold:
            print(f"SUCCESS: {test_name}")
            exit(0)
        else:
            print(f"FAIL: {test_name} at {os.getcwd()}")
            exit(1)

    if opts.verbose:

        print("max analytical {:.4e}".format(np.max(analytical_solution)))
        print("min analytical {:.4e}".format(np.min(analytical_solution)))
        print("max code {:.4e}".format(np.max(code_solution)))
        print("min code {:.4e}".format(np.min(code_solution)))

        print("integral analytical - offset {:.4e}".format(integral_analytical))
        print("integral code - offset {:.4e}".format(integral_code))
        print("integral diff {:.4e}".format(integral_diff))
        print("integral absdiff {:.4e}".format(integral_absdiff))


    # xlim=[47,53]
    # ylim=[-3,3]
    Llim = float(params["Llim"])
    xlim=ylim=[-Llim,Llim]
    # xlim = ylim = [-2, 2]
    # xlim = ylim = None
    vmax = max(np.max(analytical_solution), np.max(code_solution))

    offset = params["offset"]
    vmin = max(1e-4, offset)

    ncols = 3 if opts.solve else 2
    nrows = 2 if opts.diff else 1
    figwidth = 12 if opts.solve else 8
    figheight = 6.5 if opts.diff else 5
    fig, axsfull = plt.subplots(ncols=ncols, nrows=nrows, dpi=150,
                                sharex="all", sharey="all", figsize=(figwidth, figheight))

    if opts.diff:
        axsv = axsfull[0, :]
    else:
        axsv = axsfull

    x0 = float(params["x0"])
    Dx = g.Xi - x0
    Dy = g.Yi

    Z = analytical_solution
    # Z = Z/np.max(Z)
    fig, ax, cbar = plot_field(
        Dx, Dy, Z, scaling="log", ax=axsv[0], vmin=vmin, vmax=vmax, cmap="viridis")
    ax.set_title("analytical")
    cbar.set_label(r"$f$")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    Z = code_solution
    # Z = Z/np.max(Z)
    fig, ax, cbar = plot_field(
        Dx, Dy, Z, scaling="log", ax=axsv[1], vmin=vmin, vmax=vmax, cmap="viridis")
    ax.set_title("code")
    cbar.set_label(r"$f$")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if opts.exactmat or opts.solve:
        # check if there is a solution already saved
        # get a hash for the initial condition and DT
        h = setup_hash()
        filename = f"matexact_{h}.npy"
        if not os.path.exists(filename):
            if not opts.solve:
                print(
                    "Exact matrix solution not yet computed. Please run './check_solution -- solve' first.")
            else:
                print("Computing exact solution of linear system.")
                from solve import solve_linear_system
                with open("test_settings.yml", "r") as infile:
                    params = yaml.safe_load(infile)
                    t0 = float(params["t0"])
                initial_condition = get_solution_array("setup.yml", t0)
                Nr, Naz = initial_condition.shape
                matrix_inversion_solution = solve_linear_system(
                    outdir, initial_condition.flatten(), Nr, Naz).reshape(Nr, Naz)
                np.save(filename, matrix_inversion_solution)
        if os.path.exists(filename):
            matrix_inversion_solution = np.load(filename)
            print("max mat exact {:.4e}".format(
                np.max(matrix_inversion_solution)))
            Z = matrix_inversion_solution
            # Z = Z/np.max(Z)
            fig, ax, cbar = plot_field(
                Dx, Dy, Z, scaling="log", ax=axsv[2], vmin=vmin, vmax=vmax)
            ax.set_title("mat exact")
            cbar.set_label(r"$f$")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

    #
    # Show differences
    #

    if opts.diff:
        axsd = axsfull[1, :]

        Z = code_solution - analytical_solution
        diffmax = np.max(np.abs(Z))
        fig, ax, cbar = plot_field(
            Dx, Dy, Z, scaling="linear", ax=axsd[0], vmin=-diffmax, vmax=diffmax, cmap="bwr")
        ax.set_title("diff code - analytical")
        cbar.set_label(r"$f$")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)


    if opts.diff and opts.solve:

        Z = code_solution - matrix_inversion_solution
        diffmax = np.max(np.abs(Z))
        fig, ax, cbar = plot_field(
            Dx, Dy, Z, scaling="linear", ax=axsd[1], vmin=-diffmax, vmax=diffmax, cmap="bwr")
        ax.set_title("diff code - exact mat")
        cbar.set_label(r"$f$")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        Z = matrix_inversion_solution - analytical_solution
        diffmax = np.max(np.abs(Z))
        fig, ax, cbar = plot_field(
            Dx, Dy, Z, scaling="linear", ax=axsd[2], vmin=-diffmax, vmax=diffmax, cmap="bwr")
        ax.set_title("exact mat - analytical")
        cbar.set_label(r"$E_\mathrm{rad}$ [cgs]")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    fig.suptitle("Comparison between analytical solution and matrix solver")


    fig, axs = plt.subplots(nrows=2, height_ratios=[4, 1], dpi=150, sharex="all", figsize=(8, 6))

    fig.subplots_adjust(hspace=0.00)

    analytical_init = get_solution_array("setup.yml", float(params["t0"]))


    ax = axs[0]
    ax.plot(g.Xc[:,0]-x0, analytical_init[:,0], label="initial")
    ax.plot(g.Xc[:,0]-x0, analytical_solution[:,0], label="analytical")
    ax.plot(g.Xc[:,0]-x0, code_solution[:,0], label="code x")
    ax.set_xlim(xlim)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$f$")
    ax.legend()

    ax = axs[1]
    ax.plot(g.Xc[:,0]-x0, code_solution[:,0]/analytical_solution[:,0]-1, label="code x - analytical")
    ax.set_xlim(xlim)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\Delta f / f$")

    # ax.grid(alpha=0.5)

    if opts.outfile is not None:
        fig.savefig(opts.outfile, dpi=150, bbox_inches='tight')
    else:
        plt.show()

