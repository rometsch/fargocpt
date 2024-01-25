#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def test(out1, Nsnapshot=10, interactive=False, log=True):


    file_name1 = out1 + f"snapshots/{Nsnapshot}/Temperature1D.dat"
    data1 = np.fromfile(file_name1)
    quant1 = data1[1::4].flatten()

    fig, axs = plt.subplots(2,1,height_ratios=[3,1], sharex=True, figsize=(6,4))
    fig.subplots_adjust(hspace=0)
    
    r = data1[::4]
    r = r.flatten()

    rmax = np.max(r)
    rmin = np.min(r)
    vmin = np.min(quant1)/1.01
    vmax = np.max(quant1)*1.01

    dens = 300*np.sqrt(5/r)

    kappa = 2e-6
    nu = 5e16 # cm2/s
    sigma = 5.6704e-05 # erg cm^-2 s^-1 K^-4
    l0 = 14959787070000
    m0 = 1.98892e+33
    Sigma0 = m0 / l0 / l0
    T0 = 1.0756431684186062e+05
    G = 6.674e-8 # dyne cm^2/g^2
    omega_k = np.sqrt(G*m0*(r*l0)**(-3))
    Ttheo = np.sqrt(27/128*kappa*nu/sigma) * dens * omega_k
    Tnum = quant1 * T0

    ### Plot temperature

    ax = axs[0]
    ax2 = axs[1]
    ax.plot(r, Tnum, color="C0", label='$T$ Code', lw=2, ls="none", marker=".", markersize=5)
    ax.plot(r, Ttheo, color="k", lw=1, ls="-", label=r'$T$ Theory')

    vmin = min(vmin, np.min(Ttheo))
    vmax = max(vmax, np.max(Ttheo))

    Tdiff = np.abs(Tnum - Ttheo) / Ttheo
    axs[1].plot(r, Tdiff, label="$T$", lw=2, color="C0")


    ### Plot density

    file_name = out1 + f"snapshots/{Nsnapshot}/Sigma1D.dat"
    data_dens = np.fromfile(file_name)
    densnum = data_dens[1::4].flatten() * Sigma0

    denstheo = 300*np.sqrt(5/r)
    ax.plot(r, densnum, color="C1", label=r'$\Sigma$ Code', ls="none", marker=".", markersize=5)
    ax.plot(r, denstheo, '-k', label=r'$\Sigma$ Theory', lw=1)

    vmin2 = np.min(densnum)
    vmax2 = np.max(densnum)

    densdiff = np.abs(densnum - denstheo) / denstheo
    axs[1].plot(r, densdiff, label=r"$\Sigma$", lw=2, color="C1")
    axs[1].set_ylabel('Relative difference')

    ax.legend(loc='upper right')
    axs[1].legend(loc='lower right')

    # ax.set_xlim(rmin,rmax)
    # ax.set_ylim(vmin,vmax)
    if log:
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log")

    if log:
        ax2.set_yscale("log", nonpositive='clip')
        ax2.set_xscale("log")



    # ax.set_ylim(bottom=5)
    # ax2.set_ylim(top=700)

    ax.set_ylabel("value [cgs]")
    ax2.set_xlabel("r [au]")

    fig.savefig("plot.jpg", dpi=150, bbox_inches='tight')
    # fig.savefig("plot.pdf", dpi=300, bbox_inches='tight')

    if interactive:
        plt.show()

    # check temperature deviation
    threshold = 0.01
    rmin = 2
    rmax = 15
    radial_range = np.logical_and(r > rmin, r < rmax)
    max_diff = np.max(Tdiff[radial_range])
    pass_test = max_diff < threshold
    with open("test.log", "w") as f:
        print(f"Test name: TemperatureTest", file=f)
        print(f"Max temperature deviation: {max_diff}", file=f)
        print(f"Threshold: {threshold}", file=f)
        print(f"Radial range: {rmin} - {rmax}", file=f)
        print(f"Pass test: {pass_test}", file=f)

    if pass_test:
        print("SUCCESS: TemperatureTest")
    else:
        print("FAIL: TemperatureTest")

if __name__=="__main__":
    test("../../output/tests/TemperatureTest/out/")