#!/usr/bin/env python3

import os
import numpy as np
from types import SimpleNamespace


def main():

    test_name = "FLD1D"
    success = calc_deviation(os.getcwd()+"/output/out")

    if success:
        print(f"SUCCESS: {test_name}")
    else:
        print(f"FAIL: {test_name} at {os.getcwd()}")




def calc_deviation(outdir):

    data = load_data(outdir)

    Nfirst = data.Ns[0]
    Nlast = data.Ns[-1]


    Rmin, Rmax = 0.2, 10
    R = np.geomspace(Rmin, Rmax, 1000)
    theo = theoretical_results(data.rc)

    deltaT = data.Tprofiles[Nlast] / theo.T - 1
    dev = np.max(np.abs(deltaT[data.rc < 6]))

    success = dev < 0.1

    make_plot_temperature(data)
    make_plot_surface_density(data)
    make_plot_surface_density_2d(data)

    return success


def theoretical_results(R):
    """Calculate the theoretical equilibrium temperature profile.

    Credit: Alex Ziampras

    Args:
        R: Radius array in units of au.

    Returns:
        Namespace with results
    """
    mu = 2.353
    K = 106701.29 # code unit for temperature
    T0 = mu * 0.05**2 * K

    Rmin = R[0]
    Rmax = R[-1]
    
    f1, f2 = -3.5, 5 # old module
    f1, f2 = -2, 9/2 # new module

    R1 = Rmin ** f1
    R2 = Rmax ** f1
    T1 = (T0 / Rmin)**f2
    T2 = (T0 / Rmax)**f2
    c1 = (T2-T1) / (R2-R1)
    c2 = (R2*T1 - R1*T2) / (R2-R1)
    T = (c1 * R ** f1 + c2) ** (1/f2)

    from types import SimpleNamespace

    theo = SimpleNamespace()
    theo.R = R
    theo.T = T
    theo.T0 = T0
    theo.T0p = T0/R
    theo.Sigma = R**-0.5
    theo.H = T**0.5*R**1.5
    return theo


def load_data(outdir):

    data = SimpleNamespace()

    data.Nr, data.Naz = np.genfromtxt(
        f"{outdir}/dimensions.dat", usecols=(4, 5), unpack=True, dtype=int)
    data.ri = np.genfromtxt(f"{outdir}/used_rad.dat")  # [1:-1]
    data.phii = np.linspace(0, 2*np.pi, data.Naz+1)
    data.rc = 0.5*(data.ri[1:] + data.ri[:-1])

    data.Ns = np.genfromtxt(f"{outdir}/snapshots/list.txt", dtype=int)

    tempunit = 1
    with open(f"{outdir}/units.dat", "r") as infile:
        for line in infile:
            if len(line) > 0 and line[0] == "#":
                continue
            line = line.strip().replace("\t", " ")
            parts = line.split(" ")
            if parts[0] == "temperature":
                tempunit = float(parts[1])

    data.Ts = {n: tempunit*np.fromfile(f"{outdir}/snapshots/{n}/Temperature.dat",
                                       dtype=np.float64).reshape(data.Nr, data.Naz) for n in data.Ns}
    data.Tprofiles = {n: np.average(T, axis=1) for n, T in data.Ts.items()}

    data.Sigmas = {n: np.fromfile(f"{outdir}/snapshots/{n}/Sigma.dat",
                                  dtype=np.float64).reshape(data.Nr, data.Naz) for n in data.Ns}
    data.Sigmaprofiles = {n: np.average(S, axis=1)
                          for n, S in data.Sigmas.items()}

    return data


def make_plot_surface_density(data):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(nrows=2, dpi=150, layout="constrained")

    cmap = plt.cm.get_cmap("viridis")

    Nfirst = data.Ns[0]
    Nlast = data.Ns[-1]

    profile0 = data.Sigmaprofiles[Nfirst]
    r = 0.5*(data.ri[1:] + data.ri[:-1])

    inds = np.linspace(Nfirst, Nlast, 10, dtype=int)
    for k, n in enumerate(inds):
        color = cmap(k/(len(inds)-1))

        y = data.Sigmaprofiles[n]

        t = n

        # ax.plot(r, (profile-profile0)/profile0, label=f"t={t:.3f}yr")
        axes[0].plot(r[1:-1], y[1:-1], label=f"t={t:.0f} orb", color=color)
        y = y/profile0 - 1
        axes[1].plot(r[1:-1], y[1:-1], label=f"t={t:.0f} orb", color=color)

    axes[1].legend()
    axes[1].set_xlabel("r [au]")
    axes[0].set_ylabel(r"$\Sigma$")
    axes[1].set_ylabel(r"$\Sigma/\Sigma_0 - 1$")

    fig.savefig(os.getcwd()+"/sigma_profiles.jpg", dpi=150)


def make_plot_surface_density_2d(data):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(dpi=150, layout="constrained")

    cmap = plt.cm.get_cmap("viridis")

    Nlast = data.Ns[-1]

    Ri, Phii = np.meshgrid(data.ri, data.phii, indexing="ij")

    pcm = ax.pcolormesh(Ri, Phii, data.Sigmas[Nlast], cmap=cmap)

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(r"$\Sigma$")

    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$\phi$")

    fig.savefig(os.getcwd()+"/sigma_2d.jpg", dpi=150)


def make_plot_temperature(data):

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(nrows=2, dpi=150, layout="constrained")

    cmap = plt.cm.get_cmap("viridis")

    Nfirst = data.Ns[0]
    Nlast = data.Ns[-1]

    profile0 = data.Tprofiles[Nfirst]
    r = 0.5*(data.ri[1:] + data.ri[:-1])

    inds = np.linspace(Nfirst, Nlast, 10, dtype=int)
    for k, n in enumerate(inds):
        color = cmap(k/(len(inds)-1))

        y = data.Tprofiles[n]

        t = n

        # ax.plot(r, (profile-profile0)/profile0, label=f"t={t:.3f}yr")
        axes[0].plot(r[1:-1], y[1:-1], label=f"t={t:.0f} orb", color=color)
        y = y/profile0 - 1
        axes[1].plot(r[1:-1], y[1:-1], label=f"t={t:.0f} orb", color=color)

    axes[1].legend()
    axes[1].set_xlabel("r [au]")
    axes[0].set_ylabel(fr"$T$ [K]")
    axes[1].set_ylabel(fr"$T/T_0 - 1$")

    fig.savefig(os.getcwd()+"/temperature_profiles.jpg", dpi=150)


if __name__ == "__main__":
    main()
