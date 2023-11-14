#!/usr/bin/env python3


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("outdir", type=str, help="Output directory to load data from.")
parser.add_argument("plotfilename", type=str, help="Filename of the plot.")
opts = parser.parse_args()

import numpy as np

def calc_theo(R):
    
    mu = 2.353
    K = 106700 # code unit for temperature
    T0 = mu * 0.05**2 * K


    f1, f2 = -3.5, 5 # old module
    f1, f2 = -2, 9/2 # new module

    Rmin = R[0]
    Rmax = R[-1]

    R1 = Rmin ** f1
    R2 = Rmax ** f1
    T1 = (T0 / Rmin)**f2
    T2 = (T0 / Rmax)**f2
    c1 = (T2-T1) / (R2-R1)
    c2 = (R2*T1 - R1*T2) / (R2-R1)
    T = (c1 * R ** f1 + c2) ** (1/f2)

    from types import SimpleNamespace

    theo = SimpleNamespace()
    theo.R = R
    theo.T = T
    theo.T0 = T0
    return theo

# Calculate the radial flux to check if the divergence is close to zero.
# This is equivalent to checking that $$R \frac{T^3}{\rho} \nabla T = \text{const.}$$

def calc_flux(R, T):
    # R**1.5 = 1/rho = H/Sigma
    H = np.sqrt(T*R**3)
    sig = R ** -0.5
    rho = sig/H
    return T**3 / rho * np.gradient(T, R)

###
### Plot data
###

import matplotlib.pyplot as plt

fig, axes = plt.subplots(nrows=3, ncols=2, width_ratios=[5,1], dpi=150, figsize=(6,7), sharex="all")
empty_axes = axes[:,1]
for ax in empty_axes:
    ax.axis("off")
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
axes = axes[:,0]
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.00)



from fargocpt import Loader
l = Loader(opts.outdir)

dataunit = "J/cm2"

Nfirst = l.snapshots[0]
Nlast = l.snapshots[-1]

r, profile0 = l.gas.vars2D.avg("Temperature", Nfirst)
profile0 = profile0.to_value("K")
r = r.to_value("au")

theo = calc_theo(r)

# normalization value for flux
n = np.argmin(np.abs(theo.R - 1))
R0 = theo.R[n]
F0 = calc_flux(r, theo.T)[n]


from matplotlib import colormaps
cmap = colormaps.get_cmap("viridis")



inds = np.geomspace(np.min(1,Nfirst), Nlast, 10, dtype=int)
# inds = [0, 1, 3, 5, 8, 10, 20, 30, 50, 75, 100, Nlast]
inds = np.append([Nfirst], inds)
inds = np.unique(inds)
for k, n in enumerate(inds):
    color = cmap(k/(len(inds)-1))
    
    try:
        T = l.gas.vars2D.avg("Temperature", n, grid=False)
    except FileNotFoundError:
        continue

    y = T.to_value("K")
    
    t = l.snapshot_time[n].to_value("yr")/1**1.5
    line, = axes[0].plot(r, y, label=f"t={t:.0f} orbs", color=color)
    axes[1].plot(r, np.abs((y/theo.T-1)), label=f"t={t:.3f}yr", color=color)

    ### Plot radial flux
    F = calc_flux(r, y)
    axes[2].plot(r, r*F / (R0*F0), label=f"t={t:.3f}yr", color=color)




### Plot theoretical data

axes[0].plot(theo.R, theo.T0 / theo.R, ls=":", label='initial', color="tab:red")
axes[0].plot(theo.R, theo.T, ls=":", label='solution', color="tab:blue")

### Plot radial flux
F = calc_flux(r, theo.T)
axes[2].plot(r, r*F / (R0*F0), ls=":", label='solution', color="tab:blue")



axes[0].legend(bbox_to_anchor=(1.05,0.5))

axes[0].set_ylabel('$T$ [K]')
axes[0].set_yscale('log')

axes[1].set_ylabel('$|T/T_\mathrm{eq} - 1|$')
axes[1].set_yscale('log')

axes[2].set_ylabel(r'$\frac{r\,F}{r\,F_\mathrm{sol}}$')

axes[-1].set_xlabel('$r$ [au]')
for ax in axes:
    ax.grid(alpha=0.4)

fig.savefig(opts.plotfilename)