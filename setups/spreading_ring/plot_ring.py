"""
Test for Super-Time stepping from D'Angelo et al. 2003 THERMOHYDRODYNAMICS OF CIRCUMSTELLAR DISKS WITH HIGH-MASS PLANETS
"""
import subprocess
import os

import numpy as np
from scipy.special import iv
import matplotlib as mpl
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



def test_1D(out, ax, dt):

    Radii = np.loadtxt(out + "/used_rad.dat", skiprows=0)
    Radii = Radii
    Rinf = Radii[:-1]
    Rsup = Radii[1:]
    Rmed = 2.0/3.0*(Rsup*Rsup*Rsup-Rinf*Rinf*Rinf)
    Rmed = Rmed / (Rsup*Rsup-Rinf*Rinf)

    M = 1.0e-10
    nu = 9.841919e-06

    Quantities = np.loadtxt(out + '/Quantities.dat', skiprows=26)
    # t = Quantities[int(dt),2]
    t = Quantities[int(dt),2]
    x = Rmed
    tau0 = 0.018
    tau = 12 * nu * t + tau0

    sig0 = M / (2*np.pi) * 2 / (Radii[-2]**2 - Radii[1]**2)
    sig = 1e-4

    # print(tau, tau-tau0)

    I = iv(0.25, 2.0*x/tau)
    sig_ring = sig * I * np.exp((-1-x*x)/tau)/(np.pi * tau * pow(x, 0.25))
    Sigma = (sig_ring + sig0)
    I0 = iv(0.25, 2.0*x/tau0)
    sig_ring0 = sig * I0 * np.exp((-1-x*x)/tau0)/(np.pi * tau0 * pow(x, 0.25))
    sig_disk0 = sig0
    Sigma0 = (sig_ring0 + sig_disk0)

    ax.plot(Rmed, Sigma, ls='-', color='red', lw=2.5, label='Theory')
    ax.plot(Rmed, Sigma0, ls='--', color='black', lw=2.5, label='Initial')

    file_name = out + '/gasdens' + str(dt) + ".dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(Rmed)
    nphi = int(N/nr)

    data = data.reshape((nr, nphi))
    data_slice= np.copy(data[:, 0])

    data = np.mean(data, 1)

    # ax.vlines(1, ymin=0, ymax=1.1*np.max(Sigma))
    ax.set_title(f"tau = {tau:.3f}")
    ax.vlines(1.0, ls='--', ymin=0, ymax=1.1*np.max(Sigma0))
    ax.plot(Rmed, data_slice, ls='-', color='m', lw=1.5, label='Simulation slice')
    ax.plot(Rmed, data, ls='--', color='blue', lw=1.5, label='Simulation mean')
    # ax.axis('auto')
    # ax.set_xlim(0.6, 1.4)
    ax.legend(loc='upper right')
    return data


def test_2D(out, ax, dt):

    Radii = np.loadtxt(out + "/used_rad.dat", skiprows=0)
    Radii = Radii
    Rinf = Radii[:-1]
    Rsup = Radii[1:]
    Rmed = 2.0/3.0*(Rsup*Rsup*Rsup-Rinf*Rinf*Rinf)
    Rmed = Rmed / (Rsup*Rsup-Rinf*Rinf)

    file_name = out + '/gasdens' + str(dt) + ".dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(Rmed)
    nphi = int(N/nr)

    phi_range = np.linspace(0, 2*np.pi, nphi+1)
    phi, r = np.meshgrid(phi_range, Radii)
    X = r*np.cos(phi)
    Y = r*np.sin(phi)

    cmap = mpl.cm.get_cmap("CMRmap").copy()

    if True:
        data = data.reshape((nr, nphi))
        ax.pcolormesh(X, Y,data, cmap=cmap, linewidth=0,rasterized=True, edgecolors='black')
    else:
        data = data.reshape((nr, nphi)).T
        data = data/np.mean(data, axis=0) - 1
        ax.pcolormesh(X.T, Y.T,data, cmap=cmap, linewidth=0,rasterized=True, edgecolors='black', vmin=-0.01, vmax=0.01)

    # ax.axis('equal')
    return data

def test_3D(out, ax, dt):

    Radii = np.loadtxt(out + "/used_rad.dat", skiprows=0)
    Radii = Radii
    Rinf = Radii[:-1]
    Rsup = Radii[1:]
    Rmed = 2.0/3.0*(Rsup*Rsup*Rsup-Rinf*Rinf*Rinf)
    Rmed = Rmed / (Rsup*Rsup-Rinf*Rinf)

    file_name = out + '/gasdens' + str(dt) + ".dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(Rmed)
    nphi = int(N/nr)
    data = data.reshape((nr, nphi))

    phi_range = np.linspace(0, 2*np.pi, nphi+1)
    phi, r = np.meshgrid(phi_range, Radii)                                                                         
    X = phi
    Y = r

    cmap = mpl.cm.get_cmap("CMRmap").copy()

    ax.pcolormesh(X,Y, data, cmap=cmap)

    R0 = 1
    nu = 9.841919e-06

    Quantities = np.loadtxt(out + '/Quantities.dat', skiprows=26)
    t = Quantities[int(dt),2]
    tau0 = 0.018
    tau = 12 * nu * t / R0**2 + tau0

    time = (tau) * R0**2 / (12 * nu * 2 * np.pi)
    ax.set_title(f'Sigma, dt = {time:2f}', color='black', y = 1.06)

    return data


parfile = "log100.par"
#compile_fargo('../../')
#run('../../', 'test/spreading_ring/' + parfile)

dts = np.array([0, 12, 36, 48, 60, 84])#, 6, 14, 18, 22, 30]
dts = np.array([0, 20, 40, 60, 80, 100])#, 6, 14, 18, 22, 30]

# dts = np.array([18,54,126,162,198, 270]) - 18#, 6, 14, 18, 22, 30]
fig, axs = plt.subplots(2,3,figsize=(8,4))
fig.subplots_adjust(hspace=0.42)

axs = np.ravel(axs)


for i in range(len(dts)):
    if i > 0:
        axs[i].get_shared_x_axes().join(axs[i],axs[i-1])
        axs[i].get_shared_y_axes().join(axs[i],axs[i-1])
    dt = dts[i]
    ax = axs[i]
    data = test_1D("../../" + parfile[:-4], ax, dt)



plt.show()
