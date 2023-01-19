#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(dpi=150)

for name in ["1", "1_mpi", "2"]:
    fname = f"4cps_{name}.txt"
    nproc, nthread, ms = np.genfromtxt(fname, unpack=True)
    ncores = nproc*nthread
    ax.plot(ncores, 17/ms, label=name, marker="x")



ax.legend()

ax.set_ylabel("speedup vs 1 core")
ax.set_xlabel("Ncores")

# ax.set_xscale("log")
# ax.set_yscale("log")

ax.grid()

fig.savefig("scaling.png")