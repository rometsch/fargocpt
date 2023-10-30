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

    for i in range(len(quants)):
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

        ax.plot(r1, data, ls='-', color=color, label=label, lw=2.5)

def analytic(axs):
    analytic_data = np.loadtxt("analytic_shock.dat", skiprows=2)
    x = analytic_data[:,1]

    for i in range(len(quants)):
        ax = axs[i]
        ax.set_xlabel('x')

        if quants[i] == "pressure":
            j = 5
        if quants[i] == "Sigma":
            j = 3
        if quants[i] == "vrad":
            j = 2

        y = analytic_data[:,(j)]

        ax.plot(x, y, '-m', label='Analytic', lw=5)

def test(_):

    Nsnap = 1
    fig, axs = plt.subplots(3,1,figsize=(6,12))
    # fig, axs = plt.subplots(2,1,figsize=(6,8))
    axs = np.ravel(axs)

    analytic(axs)
    plot('../../output/tests/shocktube_PVTE/shocktube/', 'Perfect Eos FARGO', 'black', Nsnap, axs)
    plot('../../output/tests/shocktube_PVTE/shocktube_varGamm/', 'Ideal Eos FARGO', 'darkblue', Nsnap, axs)

    x0, _, rho0, v0, p0 = np.loadtxt("plutoIdealGas0.tab", unpack=True)
    x1, _, rho1, v1, p1 = np.loadtxt("plutoIdealGas1.tab", unpack=True)

    xp1, _, rhop1, vp1, pp1 = np.loadtxt("plutoPerfectGas1.tab", unpack=True)

    pluto_quants0 = np.array([rho0, v0, p0])
    pluto_quants1 = np.array([rho1, v1, p1])

    pluto_perfectQuants1 = np.array([rhop1, vp1, pp1])

    for ind in range(len(quants)):
        ax = axs[ind]

        quant = quants[ind]
        i = quants_key[quant]

        ax.axis('auto')
        ax.set_title(quant, color='black', y = 1.00)
        ax.plot(xp1, pluto_perfectQuants1[i], ls='--',color="red",  label="Perfect Eos PLUTO")
        ax.plot(x1, pluto_quants1[i],ls='--' ,color="orange", label="Ideal Eos PLUTO")
        #ax.grid()
        if quant == 'Sigma':
            ax.legend(loc='upper right')
        if quant == 'energy':
            ax.legend(loc='upper right')
        if quant == 'vrad':
            ax.legend(loc='upper left')
        if quant == 'pressure':
            ax.legend(loc='upper right')


    plt.tight_layout()
    fig.savefig('plot.jpg', dpi=150, bbox_inches='tight')
    # plt.show()

    print("NO PASSFAIL YET: shocktube_PVTE")

if __name__ == "__main__":
    test(None)