#!/usr/bin/env python3

import argparse
import numpy as np

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

Xold = np.fromfile(outdir + "/snapshots/0/energy.dat", dtype=np.float64)
a = np.fromfile(outdir + "/snapshots/1/A.dat", dtype=np.float64)
b = np.fromfile(outdir + "/snapshots/1/B.dat", dtype=np.float64)
c = np.fromfile(outdir + "/snapshots/1/C.dat", dtype=np.float64)
d = np.fromfile(outdir + "/snapshots/1/D.dat", dtype=np.float64)
e = np.fromfile(outdir + "/snapshots/1/E.dat", dtype=np.float64)

b = b.reshape(Nr, Naz)
b[0,:] = b[1,:]
b[-1,:] = b[-2,:]
b = b.flatten()

N = Nr*Naz

A = np.zeros((N,N))

for k in range(N):
    A[k,k] = b[k]
    A[k, (k-Naz)%N] = a[k]
    A[k, (k+Naz)%N] = c[k]
    A[k, (k-1)%N] = d[k]
    A[k, (k+1)%N] = e[k]
    if A[k,k] == 0:
        print(f"Diagonal element is zero for k={k}")

Xnew = np.linalg.solve(A, Xold)
print(Nr*Naz,len(Xnew))

"""
Plotting the result
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors


fig, ax = plt.subplots(dpi=150)

Z = Xnew.reshape(Nr, Naz)

print(X.shape)

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
# cbar.set_label(f"{name} [{dataunit}]")

# ax.set_xlabel(f"x [{lengthunit}]")
# ax.set_ylabel(f"y [{lengthunit}]")

plt.show()