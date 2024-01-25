import subprocess
import os

import numpy as np
import matplotlib.pyplot as plt

import yaml

quants = ["Sigma", "vrad", "pressure"]
quants_key = {'Sigma': 0, 'vrad': 1, 'pressure': 2}

def plot(out, label, color, Nsnap, axs):

    r12 = np.loadtxt(out + "used_rad.dat", skiprows=0)
    r1 = 0.5*(r12[1:] + r12[:-1])-r12[0]
    rmax = np.max(r1)
    rmin = np.min(r1)

    file_name = out + 'snapshots/' + str(Nsnap) + '/Sigma' + ".dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(r1)
    nphi = int(N/nr)

    fargo_quants = {}

    for i in range(len(quants)):
        quant = quants[i]
        ax = axs[i]

        file_name = out + 'snapshots/' + str(Nsnap) + '/' + quants[i] + ".dat"
        data = np.fromfile(file_name)
        if quants[i] == 'vrad':
            data = data.reshape((nr+1, nphi))
            data = np.mean(data, 1)
            data = 0.5*(data[1:] + data[:-1])
        else:
            data = data.reshape((nr, nphi))
            data = np.mean(data, 1)

        fargo_quants[quant] = data
    
    return r1, np.array([ fargo_quants["Sigma"], fargo_quants["vrad"], fargo_quants["pressure"]])

def test(_):

    Nsnap = 1
    fig, axs = plt.subplots(6,1,figsize=(6,12))
    # fig, axs = plt.subplots(2,1,figsize=(6,8))
    axs_diff = axs[1::2]
    axs = axs[::2]

    fargo_r, fargo_quants = plot('../../output/tests/shocktube_PVTE/shocktube_varGamm/', 'Perfect Eos FARGO', 'darkblue', Nsnap, axs)

    x1, _, rho1, v1, p1 = np.loadtxt("plutoPVTE.tab", unpack=True)

    pluto_quants = np.array([rho1, v1, p1])

    diff_quants = np.zeros_like(pluto_quants)
    for i in range(len(quants)):
        diff_quants[i] = np.abs(fargo_quants[i] - pluto_quants[i])

    for ind in range(len(quants)):
        ax = axs[ind]

        quant = quants[ind]
        i = quants_key[quant]

        ax.axis('auto')
        ax.set_title(quant, color='black', y = 1.00)

        ax.plot(fargo_r, fargo_quants[i], ls='-', color="C0", label="Perfect Eos Fargo", lw=2.5)
        ax.plot(x1, pluto_quants[i],ls='--' ,color="orange", label="Perfect Eos PLUTO")

        axs_diff[ind].plot(x1, diff_quants[i])

        #ax.grid()
        if quant == 'Sigma':
            ax.legend(loc='upper right')
        if quant == 'energy':
            ax.legend(loc='upper right')
        if quant == 'vrad':
            ax.legend(loc='upper left')
        if quant == 'pressure':
            ax.legend(loc='upper right')



    with open("testconfig.yml", "r") as f:
        testconfig = yaml.safe_load(f)
    acceptable_diffs = testconfig["acceptable_diffs"]

    success = True
    with open("test.log", "w") as logfile:
        from datetime import datetime
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{current_time}", file=logfile)

        for quant in quants:
            i = quants_key[quant]
            diff = np.average(diff_quants[i])
            threshold = acceptable_diffs[quant]

            print(f"{quant}: {diff = }, {threshold = }", file=logfile)

            is_smaller = diff < threshold
            if not is_smaller:
                success = False

    if success:
        print("SUCCESS: shocktube_PVTE")
    else:
        print("FAIL: shocktube_PVTE, threshold of averaged diff exceeded")

    plt.tight_layout()
    fig.savefig('plot.jpg', dpi=150, bbox_inches='tight')

if __name__ == "__main__":
    test(None)