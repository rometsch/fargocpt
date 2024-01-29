import subprocess
import os

import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate, integrate


import yaml

quants = ["Sigma"]
# quants = ["Sigma", "vrad", "pressure"]
quants_key = {'Sigma': 0, 'vrad': 1, 'pressure': 2}

def  plot(out, label, color, Nsnap, axs, ls=None, lw=None):

    r12 = np.loadtxt(out + "used_rad.dat", skiprows=0)
    r1 = 0.5*(r12[1:] + r12[:-1])-r12[0]
    rmax = np.max(r1)
    rmin = np.min(r1)

    file_name = out + 'snapshots/' + str(Nsnap) + '/Sigma' + ".dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(r1)
    nphi = int(N/nr)

    ax = axs[0]

    file_name = out + 'snapshots/' + str(Nsnap) + '/' + "Sigma" + ".dat"
    data = np.fromfile(file_name)
    data = data.reshape((nr, nphi))
    data = np.mean(data, 1)

    ax.plot(r1, data, color=color, label=label, lw=lw, ls=ls)
    return r1, data

def get_analytic_spl(quant):
    analytic_data = np.loadtxt("analytic_shock.dat", skiprows=2)
    x = analytic_data[:,1]
    analytic_quants_key = {"vrad": 0, "Sigma": 1, "Temperature": 2, "energy": 3}

    i = analytic_quants_key[quant]

    if i is None:
        raise ValueError(f"{quant} is not a valid quantity")
    y = analytic_data[:,(i+2)]
    if quant == 'energy':
      y = analytic_data[:,4]*analytic_data[:,3]/(1.4-1)
      
    s = interpolate.InterpolatedUnivariateSpline(x, y)

    return s

def analytic(ax, quant):
    analytic_data = np.loadtxt("analytic_shock.dat", skiprows=2)
    analytic_quants_key = {"vrad": 0, "Sigma": 1, "Temperature": 2, "energy": 3}

    x = analytic_data[:,1]
    i = analytic_quants_key[quant]
    y = analytic_data[:,(i+2)]
    if quant == 'energy':
        y = analytic_data[:,4]*analytic_data[:,3]/(1.4-1)

    ax.plot(x, y, '-k', label='Ideal analytic', lw=2)


def test(_):

    Nsnap = 1
    
    fig, axs = plt.subplots(2,1,height_ratios=[3,1], sharex=True, figsize=(6,4))
    fig.subplots_adjust(hspace=0)
    
    # fig, axs = plt.subplots(2,1,figsize=(6,8))
    axs = np.ravel(axs)

    analytic(axs[0], "Sigma")
    r_fargo_ideal, Sigma_fargo_ideal = plot('../../output/tests/shocktube_PVTE/shocktube/', 'Ideal FargoCPT', 'C0', Nsnap, axs, lw=2, ls="--")

    x1, _, rho1, v1, p1 = np.loadtxt("plutoIdeal.tab", unpack=True)

    xp1, _, rhop1, vp1, pp1 = np.loadtxt("plutoPVTE.tab", unpack=True)

    pluto_quants1 = np.array([rho1, v1, p1])

    pluto_perfectQuants1 = np.array([rhop1, vp1, pp1])

   
    ax = axs[0]

    quant = "Sigma"
    i = 0

    ax.axis('auto')
    # ax.set_title(quant, color='black', y = 1.00)
    # ax.plot(x1, pluto_quants1[i],ls='--' ,color="orange", label="Ideal Eos PLUTO")
    ax.plot(xp1, pluto_perfectQuants1[i], ls='-',color="black",  label="PVTE PLUTO")
    #ax.grid()
    
    r_fargo_pvte, Sigma_fargo_pvte = plot('../../output/tests/shocktube_PVTE/shocktube_varGamm/', 'PVTE FargoCPT', 'C1', Nsnap, axs, lw=2, ls="--")
    ax.legend(loc='upper right')
    ax.set_ylabel(r"$\Sigma$")


    ax = axs[1]
    spl = get_analytic_spl(quant)
    Sigma_ana_ideal = spl(r_fargo_ideal)
    ax.plot(r_fargo_ideal, np.abs(Sigma_fargo_ideal - Sigma_ana_ideal), color="C0", lw=2, ls="-", label="Ideal")
    
    ax.plot(r_fargo_pvte, np.abs(Sigma_fargo_pvte - rhop1), color="C1", lw=2, ls="-", label="PVTE")
    
    ax.set_yscale("log")
    ax.legend(loc="upper left")
    ax.set_ylim(bottom=1e-5)
    ax.set_ylabel(r"$\Delta \Sigma$")
    ax.set_xlabel(r"$x$")
    


    with open("testconfig.yml", "r") as f:
        testconfig = yaml.safe_load(f)
    acceptable_diffs = testconfig["acceptable_diffs"]

    # success = True
    # with open("diffs.log", "w") as logfile:
    #     from datetime import datetime
    #     current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #     print(f"{current_time}", file=logfile)

    #     for key, val in test_cases.items():
    #         outdir = val["outdir"]
    #         if not os.path.exists(outdir):
    #             continue
    #         for quant in quants:
    #             diff = diff_to_analytic(outdir, quant, 1)
    #             is_smaller = diff < acceptable_diff[quant]
    #             if not is_smaller:
    #                 success = False
    #             print("SUCCESS" if is_smaller else "FAIL", "|",
    #                   key, "|", quant, "|", diff, file=logfile)

    # plt.tight_layout()
    fig.savefig('paper-plot.jpg', dpi=150)#, bbox_inches='tight')
    fig.savefig('paper-plot.pdf', dpi=300, bbox_inches='tight')
    # plt.show()

    print("NO PASSFAIL YET: shocktube_PVTE")

if __name__ == "__main__":
    test(None)