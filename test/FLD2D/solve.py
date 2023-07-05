#!/usr/bin/env python3

import argparse
import numpy as np


"""
Solve the linear system Ax = b, give by the FLD diffusion equation.
b = energy at beginning of timestep.
A = matrix set up by diffusion eq.

Xold(i,j) = aij X(i-1,j) + cij X(i+1,j) + dij X(i,j-1) + eij X(i,j+1) + bij X(i,j)

This translates to an equation with flattened X of where k = i*Nrad + Naz
xold = A * x

with x is vector of size N = Nrad * Naz
and A is a N X N matrix

The matrix is constructed such that

A(k, k) = bij
A(k, (k-Nr)%N) = aij
A(k, (k+Nr)%N) = cij
A(k, (k-1)%N) = dij
A(k, (k+1)%N) = eij
"""


def solve_linear_system(outdir, Xold, Nr, Naz):

    a = np.fromfile(outdir + "/snapshots/1/A.dat", dtype=np.float64)
    b = np.fromfile(outdir + "/snapshots/1/B.dat", dtype=np.float64)
    c = np.fromfile(outdir + "/snapshots/1/C.dat", dtype=np.float64)
    d = np.fromfile(outdir + "/snapshots/1/D.dat", dtype=np.float64)
    e = np.fromfile(outdir + "/snapshots/1/E.dat", dtype=np.float64)

    # b = b.reshape(Nr, Naz)
    # b[0,:] = b[1,:]
    # b[-1,:] = b[-2,:]
    # b = b.flatten()

    N = Nr*Naz

    A = np.zeros((N,N))


    for k in range(N):
        if b[k] == 0:
            A[k,k] = 1
        else:
            A[k,k] = b[k]
        A[k, (k-Naz)%N] = a[k]
        A[k, (k+Naz)%N] = c[k]
        if k%Naz == 0:
            A[k, (k+Naz-1)%N] = d[k]
        else:
            A[k, (k-1)%N] = d[k]
        if (k+1)%Naz == 0:
            A[k, (k-Naz+1)%N] = e[k]
        else:
            A[k, (k+1)%N] = e[k]

    for k in range(N):
        if A[k,k] == 0:
            print(f"Zero diagonal element for k = {k}")


    Xnew = np.linalg.solve(A, Xold)

    return Xnew


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", type=str, default=None, help="Output directory")


    opts = parser.parse_args()

    outdir = opts.d if opts.d is not None else "output/out"

    Nr, Naz = np.genfromtxt(outdir + "/dimensions.dat", usecols=(4,5), dtype=int)

    ri = np.genfromtxt(outdir + "/used_rad.dat")
    phii = np.linspace(0, 2*np.pi, Naz+1)
    R, PHI = np.meshgrid(ri, phii, indexing="ij")
    X = R*np.cos(PHI)
    Y = R*np.sin(PHI)

    internal_energy = np.fromfile(outdir + "/snapshots/0/energy.dat", dtype=np.float64)
    # convert internal energy to energy density
    rho = 1 # in cgs
    cV = 8.76e7
    Trad = internal_energy/rho/cV
    aR = 1.89e-15
    Erad = aR*Trad**4

    Xold = Erad
    Xnew = solve_linear_system(outdir,Erad)

    """
    Plotting the result
    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as mplcolors

    # fig, ax = plt.subplots()
    # ax.imshow(np.log10(np.abs(A)),cmap="magma")

    # for n in range(Nr):
    #     ax.axhline(n*Naz-0.5,color="black", alpha=0.4, lw=0.1)
    #     ax.axvline(n*Naz-0.5,color="black", alpha=0.4, lw=0.1)

    # plt.show()
    # exit()

    fig, ax = plt.subplots(dpi=150)

    Z = Xnew.reshape(Nr, Naz)

    vmin = np.min(Z)
    vmax = np.max(Z)

    logplot = False
    cmap = "magma"
    if logplot:
        vmin = 1e-20 + np.min(np.abs(Z))
        vmax = 1e-20 + np.max(np.abs(Z))
        print(vmin, vmax)
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
        Z = np.abs(Z)
    else:
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)
    ax.set_title("Exact solution")

    # ax.set_xlabel(f"x [{lengthunit}]")
    # ax.set_ylabel(f"y [{lengthunit}]")


    """
    Plot exact solution symmetry
    """
    fig, ax = plt.subplots(dpi=150)

    Z = Xnew.reshape(Nr, Naz)
    Z = Z - np.roll(Z, -1, axis=1)[:,::-1]


    vmin = np.min(Z)
    vmax = np.max(Z)

    logplot = False
    cmap = "magma"
    if logplot:
        vmin = 1e-20 + np.min(np.abs(Z))
        vmax = 1e-20 + np.max(np.abs(Z))
        print(vmin, vmax)
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
        Z = np.abs(Z)
    else:
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)
    ax.set_title("Exact solution symmetry")

    """
    Plot SOR solution
    """

    fig, ax = plt.subplots(dpi=150)

    Z = np.fromfile(outdir + "/snapshots/1/Erad.dat", dtype=np.float64).reshape(Nr, Naz)

    vmin = np.min(Z)
    vmax = np.max(Z)

    logplot = False
    cmap = "magma"
    if logplot:
        vmin = 1e-20 + np.min(np.abs(Z))
        vmax = 1e-20 + np.max(np.abs(Z))
        print(vmin, vmax)
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
        Z = np.abs(Z)
    else:
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)
    ax.set_title("SOR solution")

    """
    Plot SOR solution symmetry
    """

    fig, ax = plt.subplots(dpi=150)

    Z = np.fromfile(outdir + "/snapshots/1/Erad.dat", dtype=np.float64).reshape(Nr, Naz)
    Z = Z - np.roll(Z, -1, axis=1)[:,::-1]

    vmin = np.min(Z)
    vmax = np.max(Z)

    logplot = False
    cmap = "magma"
    if logplot:
        vmin = 1e-20 + np.min(np.abs(Z))
        vmax = 1e-20 + np.max(np.abs(Z))
        print(vmin, vmax)
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
        Z = np.abs(Z)
    else:
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)
    ax.set_title("SOR solution symmetry")

    """
    Plot initial condition
    """

    fig, ax = plt.subplots(dpi=150)

    Z = Xold.reshape(Nr, Naz)

    vmin = np.min(Z)
    vmax = np.max(Z)

    logplot = False
    cmap = "magma"
    if logplot:
        vmin = 1e-20 + np.min(np.abs(Z))
        vmax = 1e-20 + np.max(np.abs(Z))
        print(vmin, vmax)
        norm = mplcolors.LogNorm(vmin=vmin, vmax=vmax)
        Z = np.abs(Z)
    else:
        norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)
    ax.set_title("initial condition")

    """
    Plot difference
    """

    fig, ax = plt.subplots(dpi=150)

    Z = np.fromfile(outdir + "/snapshots/1/Erad.dat", dtype=np.float64).reshape(Nr, Naz) - Xnew.reshape(Nr, Naz)

    vmin = np.min(Z)
    vmax = np.max(Z)

    logplot = False
    cmap = "bwr"
    norm = mplcolors.SymLogNorm(vmin=-vmax, vmax=vmax, linthresh = 1e-2*vmax)

    pcm = ax.pcolormesh(X,Y,Z, norm=norm, cmap=cmap)
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)
    ax.set_title("difference SOR - exact")


    plt.show()