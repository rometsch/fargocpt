"""
Test for Super-Time stepping from D'Angelo et al. 2003 THERMOHYDRODYNAMICS OF CIRCUMSTELLAR DISKS WITH HIGH-MASS PLANETS
"""
import subprocess
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate, integrate

quants = ["vrad", "Sigma", "Temperature", "energy"]

acceptable_diff = {
    "vrad" : 0.015, 
    "Sigma" : 0.0073,
    "Temperature" : 0.016, 
    "energy" : 0.014
}

test_cases = {
    "SN" : {
        "outdir" : 'output/SN/',
        "color" : "red",
        "ls" : "--"
    },
    "TW" : {
        "outdir" : 'output/TW/',
        "color" : "blue",
        "ls" : "--"
    },
    "SN STS" : {
        "outdir" : 'output/STS/',
        "color" : "gold",
        "ls" : ":"
    },
    "TW LF" : {
        "outdir" : 'output/TW_LF/',
        "color" : "green",
        "ls" : "-."
    },
    "SN LF" : {
        "outdir" : 'output/SN_LF/',
        "color" : "orange",
        "ls" : ":"
    }
}

def compile_fargo(fargo_path):
    wd = os.getcwd()
    os.chdir(fargo_path)

    with open("make.log", "w") as logfile:
        subprocess.run(['make', '-j', '-C' 'src/'], stdout=logfile)
    os.chdir(wd)


def run(par_file):
    logfilename = os.path.basename(par_file).replace(".yml",".log")
    with open(logfilename, "w") as logfile:
        subprocess.call('../../fargo start ' + par_file, shell=True, stdout=logfile)

def get_analytic_spl(quant):
    analytic_data = np.loadtxt("analytic_shock.dat", skiprows=2)
    x = analytic_data[:,1]
    i = None
    for n, s in enumerate(quants):
        if s == quant:
            i = n
    if i is None:
        raise ValueError(f"{quant} is not a valid quantity")
    y = analytic_data[:,(i+2)]
    if quants[i] == 'energy':
      y = analytic_data[:,4]*analytic_data[:,3]/(1.4-1)
      
    s = interpolate.InterpolatedUnivariateSpline(x, y)

    return s

def analytic(axs):
    analytic_data = np.loadtxt("analytic_shock.dat", skiprows=2)
    x = analytic_data[:,1]

    for i in range(len(quants)):
        ax = axs[i]
        y = analytic_data[:,(i+2)]
        if quants[i] == 'energy':
            y = analytic_data[:,4]*analytic_data[:,3]/(1.4-1)

        ax.plot(x, y, '-k', label='Analytic', lw=2)


def plot_output(out, label, color, Nsnap, axs, ls="--"):

    r12 = np.loadtxt(out + "used_rad.dat", skiprows=0)
    r1 = 0.5*(r12[1:] + r12[:-1])-r12[0]
    file_name = f"{out}/snapshots/{Nsnap}/Sigma.dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(r1)
    nphi = int(N/nr)

    for i in range(len(quants)):
        ax = axs[i]

        name = quants[i]

        file_name = f"{out}/snapshots/{Nsnap}/{name}.dat"
        data = np.fromfile(file_name)

        if name == 'vrad':
            data = data.reshape((nr+1, nphi))
            data = np.mean(data, axis=1)
            data = 0.5*(data[1:] + data[:-1])
        else:
            data = data.reshape((nr, nphi))
            data = np.mean(data, axis=1)


        ax.plot(r1, data, ls=ls, color=color, label=label, lw=2.5)


def diff_to_analytic(out, quant, Nsnap):

    r12 = np.loadtxt(out + "used_rad.dat", skiprows=0)
    r1 = 0.5*(r12[1:] + r12[:-1])-r12[0]
    file_name = f"{out}/snapshots/{Nsnap}/Sigma.dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(r1)
    nphi = int(N/nr)


    name = quant
    file_name = f"{out}/snapshots/{Nsnap}/{name}.dat"
    data = np.fromfile(file_name)
    if name == 'vrad':
        data = data.reshape((nr+1, nphi))
        data = np.mean(data, 1)
        data = 0.5*(data[1:] + data[:-1])
    else:
        data = data.reshape((nr, nphi))
        data = np.mean(data, 1)
    # restrict to 0 to 1
    inds = np.logical_and(r1 >= 0, r1 <=1)
    r1 = r1[inds]
    data = data[inds] 

    spl = get_analytic_spl(quant)
    analytic = spl(r1)

    delta = np.abs(data - analytic)
    return integrate.simpson(delta, r1)


def visualize(Nsnapshot):
    fig, axs = plt.subplots(2,2,figsize=(8,4))
    axs = np.ravel(axs)

    analytic(axs)
    for key, val in test_cases.items():
        if not os.path.exists(val["outdir"]):
            continue
        plot_output(val["outdir"], 
                    key, 
                    val["color"], 
                    Nsnapshot, 
                    axs, 
                    ls=val["ls"])

    for i in range(len(quants)):
        ax = axs[i]
        ax.axis('auto')
        ax.set_title(quants[i], color='black', y = 1.06)
        if quants[i] == 'dens':
            ax.legend(loc='upper right')
        if quants[i] == 'energy':
            ax.legend(loc='upper right')
        if quants[i] == 'vrad':
            ax.legend(loc='upper left')
        if quants[i] == 'Temperature':
            ax.legend(loc='upper left')
    plt.show()

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--viz", action="store_true", help="Plot the results.")
    opt = parser.parse_args()
    
    compile_fargo('../../')
    run('setups/shocktube_SN.yml')
    run('setups/shocktube_TW.yml')
    run('setups/shocktube_STS.yml')
    run('setups/shocktube_TW_LF.yml')
    run('setups/shocktube_SN_LF.yml')
    
    success = True
    with open("diffs.log", "w") as logfile:
        for key, val in test_cases.items():
            outdir = val["outdir"]
            if not os.path.exists(outdir):
                continue
            for quant in quants:
                diff = diff_to_analytic(outdir, quant, 228)
                is_smaller = diff < acceptable_diff[quant]
                if not is_smaller:
                    success = False
                print("SUCCESS" if is_smaller else "FAIL", "|",
                      key, "|", quant, "|", diff, file=logfile)
    if success:
        print("SUCCESS: shocktube test")
    else:
        print("FAIL: shocktube test, threshold of integrated diff exceeded")
    
    if opt.viz:
       visualize(228)

if __name__=="__main__":
    main()