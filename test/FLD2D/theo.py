#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

# below, all physical quantities are in cgs

energy0 = 1
K = 1e10 # diffusion coefficient

def energy(t, x, y):

    r = np.sqrt(x**2 + y**2)

    return energy0 / (4*np.pi*K*t) * np.exp( - r**2/(4*K*t))


L = -10
Ncells = 100
xi = np.linspace(-L, L, Ncells+1)
yi = np.linspace(-L, L, Ncells+1)
x = 0.5*(xi[1:] + xi[:-1])
y = 0.5*(yi[1:] + yi[:-1])
Xi, Yi = np.meshgrid(xi, yi, indexing="ij")
X, Y = np.meshgrid(x, y, indexing="ij")

t = 1e-10

fig, ax = plt.subplots()

pcm = ax.pcolormesh(Xi, Yi, energy(t, X, Y), cmap="magma")

cbar = fig.colorbar(pcm, ax=ax)
ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
cbar.set_label("energy [erg]")
ax.set_title(f"t = {t:.3e} sec")

plt.show()