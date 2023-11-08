#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
import astropy.units as u

# below, all physical quantities are in cgs

energy0 = 1*u.erg
mu = 2.35
gamma = 1.4
rho = 1 * u.g/u.cm**3
kappa = 1 * u.cm**2/u.g

Rg = const.k_B / const.m_p
cv = Rg / mu / (gamma - 1)
cv = cv.cgs

sigmaR = const.sigma_sb
aR = 4*sigmaR/const.c
aR = aR.to("erg/(cm3*K4)")

K = 1/3 * const.c / rho / kappa
K = K.cgs

print(f"K = {K:.3e}")

def radiative_energy(t, x, y):

    r = np.sqrt(x**2 + y**2)

    return energy0 / (4*np.pi*K*t) * np.exp( - r**2/(4*K*t))


Sigma = 1 * u.g/u.cm**2
H = 1 * u.cm

def internal_energy_density(erad):
    Trad = (erad/H / aR)**(1/4)
    print(Trad)
    Trad = Trad.to("K")
    Tgas = Trad
    egas = Sigma*cv*Tgas
    return egas.to("erg/cm2")

L = 10*u.cm
Ncells = 100
xi = np.linspace(-L, L, Ncells+1)
yi = np.linspace(-L, L, Ncells+1)
x = 0.5*(xi[1:] + xi[:-1])
y = 0.5*(yi[1:] + yi[:-1])
Xi, Yi = np.meshgrid(xi, yi, indexing="ij")
X, Y = np.meshgrid(x, y, indexing="ij")

t = 1e-10*u.s

fig, ax = plt.subplots()

erad = radiative_energy(t, X, Y).to("erg/cm2")
print(erad/H)
egas = internal_energy_density(erad).to_value("erg/cm2")

pcm = ax.pcolormesh(Xi.to_value("cm"), Yi.to_value("cm"), egas, cmap="magma")

cbar = fig.colorbar(pcm, ax=ax)
ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
cbar.set_label("energy [erg]")
ax.set_title(f"t = {t:.3e} sec")

plt.show()