#!/usr/bin/env python3
# nodes \t ms
import numpy as np
import matplotlib.pyplot as plt

nodes = np.array([1,2,4,8,16])
cores = nodes*28
time = np.array([130,65,30,16.5,12])

fig, ax = plt.subplots()

ax.plot(cores, 1 / (time/time[0]*cores/cores[0]), marker="*", markersize=10)

ax.set_xlabel("number of cores")
ax.set_ylabel("speedup at constant work / strong scaling")

ax.set_xscale("log")
# ax.set_yscale("log")

ax.set_xticks([])
ax.set_xticks(cores)
ax.set_xticklabels([f"{c}" for c in cores])

# yticks = [10,20,30,40,50,60,70,80,90,100,110,120,130]
# ax.set_yticks(yticks)
# ax.set_yticklabels([f"{y}" for y in yticks])

ax.grid()

fig.savefig("scaling.svg")

plt.show()