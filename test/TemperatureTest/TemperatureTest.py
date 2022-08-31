#!/usr/bin/env python3
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

    subprocess.call('mpirun -n 4 ./fargo start ' + par_file, shell=True)
    os.chdir(wd)

def test(out1, dt):


    file_name1 = out1 + 'gasTemperature1D' + str(dt) + ".dat"
    data1 = np.fromfile(file_name1)
    print(file_name1)
    r1 = data1[::4]
    quant1 = data1[1::4]

    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    rmax = np.max(r1)
    rmin = np.min(r1)
    vmin = np.min(quant1)/1.01
    vmax = np.max(quant1)*1.01


    r = data1[::4]
    r = r.flatten()
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
    ax.axis('auto')
    ax.set_title('Temperature', color='black', y = 1.06)
    T = np.sqrt(27/128*kappa*nu/sigma) * dens * omega_k
    ax.plot(r, T, '--k', label='Theory', lw=2.5)

    vmin = min(vmin, np.min(T))
    vmax = max(vmax, np.max(T))

    ax.plot(r1.flatten(), quant1.flatten() * T0, '.r', label='Code', lw=2.5)

    file_name = out1 + 'gasdens1D' + str(dt) + ".dat"
    data_dens = np.fromfile(file_name)
    quant2 = data_dens[1::4] * Sigma0

    ax2.axis('auto')
    ax2.set_title('Density', color='black', y = 1.06)
    dens = 300*np.sqrt(5/r1.flatten())
    ax2.plot(r1.flatten(), dens, '--k', label='Theory', lw=2.5)
    ax2.plot(r1.flatten(), quant2, '.r', label='Code', lw=2.5)

    vmin2 = np.min(quant2)
    vmax2 = np.max(quant2)


    # ax.plot(r2.flatten(), quant2.flatten(), '--b', label='Expl', lw=2)



    ax.legend(loc='upper right')
    ax2.legend(loc='upper right')


    log = True
    ax.set_xlim(rmin,rmax)
    ax.set_ylim(vmin,vmax)
    if log:
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log")

    ax2.set_xlim(rmin,rmax)
    ax2.set_ylim(vmin2,vmax2)
    if log:
        ax2.set_yscale("log", nonpositive='clip')
        ax2.set_xscale("log")



    plt.show()

compile_fargo('../../')
run('../../', 'test/TemperatureTest/angelo.par')
test('../../angelo/', 10)
