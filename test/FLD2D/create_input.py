#!/usr/bin/env python3

import numpy as np
import yaml
import astropy.constants as const
import astropy.units as units
from types import SimpleNamespace

def get_fargo_grid(setupfile):

    # get size of grid specified in the setup file
    with open(setupfile, "r") as infile:
        params = yaml.safe_load(infile)
        Nrad = params["Nrad"]
        Naz = params["Nsec"]
        rmin = params["Rmin"]
        rmax = params["Rmax"]

    ri = np.linspace(rmin, rmax, Nrad+1)
    phii = np.linspace(0, 2*np.pi, Naz+1)
    Ri, Phii = np.meshgrid(ri, phii, indexing="ij")
    Xi = Ri*np.cos(Phii)
    Yi = Ri*np.sin(Phii)

    rc = 2/3*(ri[1:]**2/(ri[1:]+ri[:-1])) + ri[:-1] # approx center in polar coords
    phic = 0.5*(phii[1:]+phii[:-1])
    Rc, Phic = np.meshgrid(rc, phic, indexing="ij")
    Xc = Rc*np.cos(Phic)
    Yc = Rc*np.sin(Phic)

    dphi = phii[1] - phii[0]
    dr = ri[1:] - ri[:-1]
    A = 0.5*(Ri[1:,1:]**2 - Ri[:-1,1:]**2)*dphi

    return SimpleNamespace(
        Nrad=Nrad, Naz=Naz,
        ri=ri, phii=phii, Ri=Ri, Phii=Phii, Xi=Xi, Yi=Yi,
        rc=rc, phic=phic, Rc=Rc, Phic=Phic, Xc=Xc, Yc=Yc, 
        dphi=dphi, dr=dr, A=A
    )

def analytical_solution(r, t, E0, K, offset=0):
    return 3*E0/(4*np.pi*t*K)*np.exp(-r**2/(4*K*t)) + offset

def get_setup_params():
    with open("test_settings.yml", "r") as infile:
        params = yaml.safe_load(infile)
    return params

def get_solution_array(setupfile, t, offset=0):
    g = get_fargo_grid(setupfile)

    nr = g.Nrad//2
    nphi = 0

    # print(dr[nr])

    c = 2.997e10
    rho = 1
    kappa = 1
    lam = 1/3
    K = lam*c/(rho*kappa)

    print(f"Diffusion coefficient K = {K} = {K:e}")

    params = get_setup_params()

    E0 = float(params["E0"])
    energy0 = E0 / g.A[nr, nphi]
    offset = float(params["offset"])*energy0

    xcell = g.Xc[nr, nphi]
    ycell = g.Yc[nr, nphi]
    DX = g.Xc - xcell
    DY = g.Yc - ycell
    Dist = np.sqrt(DX**2 + DY**2)

    return analytical_solution(Dist,t, energy0, K, offset=offset)

def Erad_to_Eint(Erad):
    """ Convert radiative energy density to internal energy density based on the test model. """
    mu = 2.35
    gamma = 1.4

    Rg = const.k_B / const.m_p
    cv = Rg / mu / (gamma - 1)
    cv = cv.cgs.value
    Sigma = 1 # g/cm3

    sigmaR = const.sigma_sb
    aR = 4*sigmaR/const.c
    aR = aR.to_value("erg/(cm3*K4)")
    Trad = (Erad/aR)**0.25

    # assume instantaneous equilibrium
    Tgas = Trad

    Eint = cv*Sigma*Tgas
    return Eint



if __name__ == "__main__":
    # now create an energy array with values in cgs with the energy located in one cell
    
    # get initial time
    with open("test_settings.yml", "r") as infile:
        params = yaml.safe_load(infile)
        t0 = float(params["t0"])

    Erad = get_solution_array("setup.yml", t0, offset=1e5)
    initial_condition = Erad_to_Eint(Erad)

    print("Erad size", Erad.shape)
    print("initial_condition size", initial_condition.shape)

    np.array(Erad, dtype=np.float64).tofile("output/out/Erad_input.dat")
    # save this energy array to file
    initial_condition.tofile("output/out/snapshots/0/energy.dat")

    import os
    if os.path.exists("output/out/snapshots/damping"):
        initial_condition.tofile("output/out/snapshots/damping/energy.dat")

