"""
Test for Super-Time stepping from D'Angelo et al. 2003 THERMOHYDRODYNAMICS OF CIRCUMSTELLAR DISKS WITH HIGH-MASS PLANETS
"""
import subprocess
import os

import numpy as np
import matplotlib.pyplot as plt


def compile_fargo(fargo_path):
    wd = os.getcwd()
    os.chdir(fargo_path)

    subprocess.run(['make', '-j', '-C' 'src/'])
    os.chdir(wd)


def run(fargo_path, par_file):
    wd = os.getcwd()
    os.chdir(fargo_path)

    subprocess.call('./fargo start ' + par_file, shell=True)
    os.chdir(wd)


quants = ["vrad", "Sigma", "Temperature", "energy"]

def analytic(axs):
    analytic_data = np.loadtxt("analytic_shock.dat", skiprows=2)
    x = analytic_data[:,1]

    for i in range(len(quants)):
        ax = axs[i]
        y = analytic_data[:,(i+2)]
        if quants[i] == 'energy':
            y = analytic_data[:,4]*analytic_data[:,3]/(1.4-1)

        ax.plot(x, y, '-k', label='Analytic', lw=2)


def test(out, label, color, dt):

    r12 = np.loadtxt(out + "used_rad.dat", skiprows=0)
    r1 = 0.5*(r12[1:] + r12[:-1])-r12[0]
    rmax = np.max(r1)
    rmin = np.min(r1)

    file_name = out + 'snapshots/' + str(dt) + '/Sigma' + ".dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(r1)
    nphi = int(N/nr)

    for i in range(len(quants)):
        ax = axs[i]

        file_name = out + '/snapshots/' + str(dt) + '/' + quants[i] + ".dat"
        data = np.fromfile(file_name)
        if quants[i] == 'vrad':
            data = data.reshape((nr+1, nphi))
            data = np.mean(data, 1)
            data = 0.5*(data[1:] + data[:-1])
        else:
            data = data.reshape((nr, nphi))
            data = np.mean(data, 1)

        ax.plot(r1, data, ls='--', color=color, label=label, lw=2.5)




compile_fargo('../../')
run('../../', 'test/shockTube/shocktube_TW.yml')
run('../../', 'test/shockTube/shocktube_SN.yml')

dt = 228
fig, axs = plt.subplots(2,2,figsize=(8,4))
axs = np.ravel(axs)

analytic(axs)
test('../../shocktube_TW/', 'TW', 'red', dt)
test('../../shocktube_SN/', 'SN', 'blue', dt)

for i in range(len(quants)):
    ax = axs[i]
    ax.axis('auto')
    ax.set_title(quants[i], color='black', y = 1.06)
    if quants[i] == 'Sigma':
        ax.legend(loc='upper right')
    if quants[i] == 'energy':
        ax.legend(loc='upper right')
    if quants[i] == 'vrad':
        ax.legend(loc='upper left')
    if quants[i] == 'Temperature':
        ax.legend(loc='upper left')



plt.show()
