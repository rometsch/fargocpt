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

    ri = np.genfromtxt("output/out/used_rad.dat")
    phii = np.linspace(0, 2*np.pi, Naz+1)
    Ri, Phii = np.meshgrid(ri, phii, indexing="ij")
    Xi = Ri*np.cos(Phii)
    Yi = Ri*np.sin(Phii)

    rc = 2/3*(ri[1:]**2/(ri[1:]+ri[:-1]) + ri[:-1]) # approx center in polar coords
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

def analytical_solution(r, t, x0, K, offset=0.0):
    return x0/(4*np.pi*t*K)*np.exp(-(r**2)/(4*K*t)) + offset

def get_setup_params():
    with open("test_settings.yml", "r") as infile:
        params = yaml.safe_load(infile)
    return params

def get_solution_array(setupfile, t):
    g = get_fargo_grid(setupfile)
    params = get_setup_params()

    f0 = float(params["f0"])
    finit = f0 #/ g.A[nr, nphi]
    offset = float(params["offset"])*finit


    # determine location of cell closest to x0
    x0 = params["x0"]
    nr = np.argmin(np.abs(g.rc - x0))
    # nr = len(g.rc)//2
    nphi = 0

    # print(dr[nr])

    K = params["K"]
    # print(f"Analytical solution at t = {t:e}, K = {K}, y0 = {finit}, offset = {offset}")
    
    xcell = g.Xc[nr, nphi]
    ycell = g.Yc[nr, nphi]

    DX = g.Xc - xcell
    DY = g.Yc - ycell
    Dist = np.sqrt(DX**2 + DY**2)

    return analytical_solution(Dist, t, finit, K, offset=offset)



if __name__ == "__main__":
    # now create an energy array with values in cgs with the energy located in one cell
    
    # get initial time
    with open("test_settings.yml", "r") as infile:
        params = yaml.safe_load(infile)
        t0 = float(params["t0"])

    f0 = get_solution_array("setup.yml", t0)

    print("f grid size", f0.shape)
    print("f max", np.max(f0))

    np.array(f0, dtype=np.float64).tofile("output/out/f_FLD2Dtest_input.dat")

    
