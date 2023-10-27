#!/usr/bin/env python3
# %%
import subprocess
import os

import numpy as np
import matplotlib.pyplot as plt


def compile_fargo(fargo_path):
    wd = os.getcwd()
    os.chdir(fargo_path)

    subprocess.run(['make', '-j 2', '-C' 'src/'])
    os.chdir(wd)


def run(fargo_path, par_file):
    wd = os.getcwd()
    os.chdir(fargo_path)

    subprocess.call('./run_fargo start ' + par_file, shell=True)
    os.chdir(wd)

def test(out1, dt):


    file_name1 = out1 + f"snapshots/{dt}/Temperature1D.dat"
    data1 = np.fromfile(file_name1)
    print(file_name1)
    r1 = data1[::4]
    quant1 = data1[1::4]

    fig, axs = plt.subplots(1,1, figsize=(6,4.5))
    ax2 = axs
    ax = ax2.twinx()

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
    # ax.set_title('Temperature', color='black', y = 1.06)
    T = np.sqrt(27/128*kappa*nu/sigma) * dens * omega_k
    T = 104 * (5/r1.flatten())**2
    ax.plot(r, T, '-k', label='T theory', lw=5)
    ax.set_ylabel('log10 T [K]')

    vmin = min(vmin, np.min(T))
    vmax = max(vmax, np.max(T))

    ax.plot(r1.flatten(), quant1.flatten() * T0, '.', color='C0', label='T code', lw=2.5)

    file_name = out1 + f"snapshots/{dt}/Sigma1D.dat"
    data_dens = np.fromfile(file_name)
    quant2 = data_dens[1::4] * Sigma0

    ax2.axis('auto')
    # ax2.set_title('Density', color='black', y = 1.06)
    dens = 300*np.sqrt(5/r1.flatten())
    ax2.plot(r1.flatten(), dens, '-k', label='$\Sigma$ theory', lw=5)
    ax2.plot(r1.flatten(), quant2, '.', color='C1', label='$\Sigma$ code', lw=2.5)
    ax2.set_ylabel('log10 $\Sigma$ [g/cm²]')
    ax.set_xlabel('log10 R [au]')

    vmin2 = np.min(quant2)
    vmax2 = np.max(quant2)


    # ax.plot(r2.flatten(), quant2.flatten(), '--b', label='Expl', lw=2)



    ax2.legend(loc='lower left')
    ax.legend(loc='upper right')


    log = True
    ax.set_xlim(rmin*0.95,rmax*1.05)
    ax.set_ylim(vmin*800,vmax*1.2)
    if log:
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log")

    ax2.set_ylim(vmin2*0.95,vmax2*1.2)
    if log:
        ax2.set_yscale("log", nonpositive='clip')
        ax2.set_xscale("log")



    plt.savefig('TemperatureTest.pdf', dpi=300, bbox_inches='tight')
    plt.show()

#compile_fargo('../../')
#run('../../', 'test/TemperatureTest/angelo.yml')
test('../../angelo/', 10)
