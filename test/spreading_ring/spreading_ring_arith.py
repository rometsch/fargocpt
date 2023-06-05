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

    subprocess.call('./run_fargo start ' + par_file, shell=True)
    os.chdir(wd)



def test_1D(out, ax, dt):

    Radii = np.loadtxt(out + "used_rad.dat", skiprows=0)
    Rinf = Radii[:-1]
    Rsup = Radii[1:]
    Rmed = 0.5*(Rinf + Rsup)


    R0 = 1
    R0_id = np.argmin(np.abs(Rmed-R0))
    M = 1.0
    nu = 4.77e-5 / R0**2

    Quantities = np.loadtxt(out + '/monitor/Quantities.dat', skiprows=26)
    t = Quantities[int(dt*10),2]
    x = Rmed / R0
    tau0 = 0.016
    tau = 12 * nu * t / R0**2 + tau0

    # print(tau, tau-tau0)

    I = iv(0.25, 2.0*x/tau)
    Sigma = M / (np.pi * R0**2) / tau / x**(1/4) * I * np.exp(-(1+x**2)/tau)
    I0 = iv(0.25, 2.0*x/tau0)
    Sigma0 = M / (np.pi * R0**2) / tau0 / x**(1/4) * I0 * np.exp(-(1+x**2)/tau0)

    ax.plot(Rmed, Sigma, ls='-', color='red', lw=2.5, label='Theory')
    ax.plot(Rmed, Sigma0, ls='--', color='black', lw=2.5, label='Initial')

    file_name = out + f'/snapshots/{dt}/Sigma.dat'
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(Rmed)
    nphi = int(N/nr)

    data = data.reshape((nr, nphi))
    data_slice= np.copy(data[:, 0])

    data = np.mean(data, 1)


    # ax.vlines(1, ymin=0, ymax=1.1*np.max(Sigma))
    ax.vlines(Rmed[R0_id], ls='--', ymin=0, ymax=1.1*np.max(Sigma0))
    ax.plot(Rmed, data_slice, ls='-', color='m', lw=1.5, label='Simulation slice')
    ax.plot(Rmed, data, ls='--', color='blue', lw=1.5, label='Simulation mean')
    # ax.axis('auto')
    # ax.set_xlim(0.6, 1.4)
    ax.legend(loc='upper right')
    return data


def test_2D(out, ax, dt):

    Radii = np.loadtxt(out + "used_rad.dat", skiprows=0)
    Rinf = Radii[:-1]
    Rsup = Radii[1:]
    Rmed = 0.5*(Rinf + Rsup)

    file_name = out + f'/snapshots/{dt}/Sigma.dat'
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

    Radii = np.loadtxt(out + "used_rad.dat", skiprows=0)
    Rinf = Radii[:-1]
    Rsup = Radii[1:]
    Rmed = 0.5*(Rinf + Rsup)

    file_name = out + f'/snapshots/{dt}/Sigma.dat'
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
    nu = 4.77e-5 / R0**2

    Quantities = np.loadtxt(out + '/monitor/Quantities.dat', skiprows=26)
    t = Quantities[int(dt*10),2]
    tau0 = 0.016
    tau = 12 * nu * t / R0**2 + tau0

    time = (tau) * R0**2 / (12 * nu * 2 * np.pi)
    ax.set_title(f'Sigma, t = {time:.2f}', color='black', y = 1.06)

    return data


parfile = "speith2003.yml"
compile_fargo('../../')
run('../../', 'test/spreading_ring/' + parfile)

dts = np.array([7, 17, 30, 43, 70, 86])*2-7

fig, axs = plt.subplots(2,3,figsize=(8,4))
fig.subplots_adjust(hspace=0.42)

axs = np.ravel(axs)


for i in range(len(dts)):
    if i > 0:
        axs[i].get_shared_x_axes().join(axs[i],axs[i-1])
        axs[i].get_shared_y_axes().join(axs[i],axs[i-1])
    dt = dts[i]
    ax = axs[i]
    data = test_2D_R_PHI('../../speith2003/', ax, dt)



plt.show()
