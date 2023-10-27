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

    mu = 2.35
    m_H = 1.66054e-24
    k_B = 1.38065e-16
    gamma = 1.4
    l0 = 14959787070000
    rcgs = r * l0
    m0 = 1.98847e+33
    T0 = 106700.1843026118
    G = 6.6743e-08 # dyne cm^2/g^2

    eta = 2/7
    eps = 0.5
    Rs = 4.6505e-05 * l0

    Ts = 100000

    # note missing factor gamma from eq.(16) fron D'Angelo & Marzari 2012
    # difference comes from them using the adiabatic sound speed for their scaleheight and our code
    # uses the locally isothermal soundspeed for h
    h = (eta * (1 - eps) * (k_B * Ts / (mu * m_H))**4 * (Rs / (G * m0))**4 * (rcgs / Rs)**2)**(1/7)
    WG = 0.4 * (Rs / rcgs) + h * eta

    T = Ts * np.sqrt(Rs /rcgs) * ((1-eps)*WG)**(1/4)

    ax.axis('auto')
    # ax.set_title('Temperature', color='black', y = 1.06)

    ax.plot(r, T, '-b', label='T theory', lw=5)
    ax.set_ylabel('log10 T [K]')

    vmin = min(vmin, np.min(T))
    vmax = max(vmax, np.max(T))

    ax.plot(r1.flatten(), quant1.flatten() * T0, '.', color='C0', label='T code', lw=2.5)

    file_name = out1 + f"snapshots/{dt}/aspectratio1D.dat"
    data_dens = np.fromfile(file_name)
    quant2 = data_dens[1::4]

    print(quant2 / h)

    ax2.axis('auto')
    # ax2.set_title('Density', color='black', y = 1.06)
    ax2.plot(r, h, '-r', label='$h$ theory', lw=5)
    ax2.plot(r, quant2, '.', color='C1', label='$h$', lw=2.5)
    ax2.set_ylabel('h')
    ax.set_xlabel('log10 R [au]')

    vmin2 = np.min(quant2)
    vmax2 = np.max(quant2)


    # ax.plot(r2.flatten(), quant2.flatten(), '--b', label='Expl', lw=2)



    ax2.legend(loc='center left')
    ax.legend(loc='center right')


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



    #plt.savefig('TemperatureTest.pdf', dpi=300, bbox_inches='tight')
    plt.show()

compile_fargo('../../')
run('../../', 'test/TemperatureTest/angelo.yml')
test('../../angelo/', 10)
