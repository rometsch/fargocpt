#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def test(out1, Nsnapshot=10, interactive=False, log=True):


    file_name1 = out1 + f"snapshots/{Nsnapshot}/Temperature1D.dat"
    data1 = np.fromfile(file_name1)
    quant1 = data1[1::4].flatten()

    fig, axs = plt.subplots(2,2, gridspec_kw={'height_ratios': [1, 1]})
    ax = axs[0,0]
    ax1diff = axs[1,0]
    ax2 = axs[0,1]
    ax2diff = axs[1,1]

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

    ax.axis('auto')
    ax.set_title('Temperature', color='black', y = 1.06)
    ax.plot(r, Tnum, '.r', label='Code', lw=2.5)
    ax.plot(r, Ttheo, '--k', label='Theory', lw=2.5)

    vmin = min(vmin, np.min(Ttheo))
    vmax = max(vmax, np.max(Ttheo))

    Tdiff = np.abs(Tnum - Ttheo) / Ttheo
    ax1diff.plot(r, Tdiff)
    ax1diff.set_ylabel('Relative difference')

    ### Plot density

    file_name = out1 + f"snapshots/{Nsnapshot}/Sigma1D.dat"
    data_dens = np.fromfile(file_name)
    densnum = data_dens[1::4].flatten() * Sigma0

    ax2.axis('auto')
    ax2.set_title('Density', color='black', y = 1.06)
    denstheo = 300*np.sqrt(5/r)
    ax2.plot(r, densnum, '.r', label='Code', lw=2.5)
    ax2.plot(r, denstheo, '--k', label='Theory', lw=2.5)

    vmin2 = np.min(densnum)
    vmax2 = np.max(densnum)

    densdiff = np.abs(densnum - denstheo) / denstheo
    ax2diff.plot(r, densdiff)
    ax2diff.set_ylabel('Relative difference')

    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')

    # ax.set_xlim(rmin,rmax)
    # ax.set_ylim(vmin,vmax)
    if log:
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log")

    if log:
        ax1diff.set_yscale("log", nonpositive='clip')
        ax1diff.set_xscale("log")

    # ax2.set_xlim(rmin,rmax)
    # ax2.set_ylim(vmin2,vmax2)
    if log:
        ax2.set_yscale("log", nonpositive='clip')
        ax2.set_xscale("log")

    if log:
        ax2diff.set_yscale("log", nonpositive='clip')
        ax2diff.set_xscale("log")

    # ax.set_ylim(bottom=5)
    # ax2.set_ylim(top=700)

    fig.savefig("plot.jpg", dpi=150)
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
