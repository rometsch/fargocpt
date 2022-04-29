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


quants = ["gasdens", "gasvrad", "gaspressure"]

x0, _, rho0, v0, p0 = np.loadtxt("data.0000.tab", unpack=True)
x1, _, rho1, v1, p1 = np.loadtxt("data.0001.tab", unpack=True)

xp1, _, rhop1, vp1, pp1 = np.loadtxt("plutoPerfectGas1.tab", unpack=True)

pluto_quants0 = np.array([rho0, v0, p0])
pluto_quants1 = np.array([rho1, v1, p1])

pluto_perfectQuants1 = np.array([rhop1, vp1, pp1])

def test(out, label, color, dt):

    r12 = np.loadtxt(out + "used_rad.dat", skiprows=0)
    r1 = 0.5*(r12[1:] + r12[:-1])-r12[0]
    rmax = np.max(r1)
    rmin = np.min(r1)

    file_name = out + '/gasdens' + str(dt) + ".dat"
    data = np.fromfile(file_name)
    N = len(data)
    nr = len(r1)
    nphi = int(N/nr)

    for i in range(len(quants)):
        ax = axs[i]

        file_name = out + quants[i] + str(dt) + ".dat"
        data = np.fromfile(file_name)
        if quants[i] == 'gasvrad':
            data = data.reshape((nr+1, nphi))
            data = np.mean(data, 1)
            data = 0.5*(data[1:] + data[:-1])
        else:
            data = data.reshape((nr, nphi))
            data = np.mean(data, 1)

        ax.plot(r1, data, ls='-', color=color, label=label, lw=2.5)



compile_fargo('../../')
run('../../', 'test/shockTube/shocktube.par')
run('../../', 'test/shockTube/shocktube_varGamm.par')

dt = 200#228
fig, axs = plt.subplots(1,3,figsize=(18,6))
axs = np.ravel(axs)


test('../../shocktube/', 'Perfect Eos FARGO', 'black', dt)
test('../../shocktube_varGamm/', 'Ideal Eos FARGO', 'darkblue', dt)


for i in range(len(quants)):
    ax = axs[i]
    ax.axis('auto')
    ax.set_title(quants[i], color='black', y = 1.06)
    ax.plot(xp1, pluto_perfectQuants1[i], ls='--',color="red",  label="Perfect Eos PLUTO")
    ax.plot(x1, pluto_quants1[i],ls='--' ,color="orange", label="Ideal Eos PLUTO")
    ax.grid()
    if quants[i] == 'gasdens':
        ax.legend(loc='upper right')
    if quants[i] == 'gasenergy':
        ax.legend(loc='upper right')
    if quants[i] == 'gasvrad':
        ax.legend(loc='upper left')
    if quants[i] == 'gasTemperature':
        ax.legend(loc='upper left')


plt.tight_layout()
plt.savefig("/home/wehner/DisksInCloseBinaries/plots/adiabaticShocktube/EosComparison.png", dpi=500)

plt.show()
