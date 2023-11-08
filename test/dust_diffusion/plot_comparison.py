import h5py
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

def load_dset(filename, dset, N=None):
    with h5py.File(filename, "r") as f:
        if N is not None:
            data = f[dset][N,:]
        else:
            data = f[dset][:]
        try:
            unit = f[dset].attrs["unit"]
            data = data*u.Unit(unit)
        except KeyError:
            pass
    return data

from load_dust import get_time, get_sigma_dust
datadir = "output/dust_diffusion/"

try:
    times = load_dset("disklab.hdf5", "time").to_value("yr")
except FileNotFoundError:
    from disklab_model import main as disklab_main
    disklab_main()
    times = load_dset("disklab.hdf5", "time").to_value("yr")
rs = load_dset("disklab.hdf5", "r").to_value("au")


fig, ax = plt.subplots(dpi=150)

toffset = 0

# inds = range(1,100,20)
# inds = [1,5,10,12]
inds = [0, 1]
Ymax = 0
n = 1
t = get_time(datadir, n)
sigma_dust, rmid, dr = get_sigma_dust(datadir, n, nbins=101)
factor = np.sum(sigma_dust*rmid*dr*2*np.pi)
Y = sigma_dust/factor
line = ax.plot(rmid, Y, drawstyle="steps-mid", label=f"t = {t:.3f}")


Ymax = max(Ymax, np.max(Y))

n_ana = np.argmin(np.abs(times-t))
sigmad = load_dset("disklab.hdf5", "sigmad", N=n_ana)
print(sigmad.shape)
print(rs.shape)

inds = np.logical_and(rs > 7, rs < 15)
data = np.stack([rs[inds], sigmad[inds].to_value("g/cm2")], axis=1)
np.savetxt("reference_data.txt", data)

sigmad_ref = np.interp(rmid, rs, sigmad.to_value("g/cm2"))

dr = np.pad(rmid[1:]-rmid[:-1], (0,1), mode="edge")
print(np.sum(rmid*dr*2*np.pi))
factor = np.sum(sigmad_ref*rmid*dr*2*np.pi)
Y = sigmad_ref/factor




ax.plot(rmid, Y, color="C1", ls="-", label=f"analytic t = {times[n_ana]:3g} yr")
Ymax = max(Ymax, np.max(Y))
    
ax.set_yscale("log")
ax.set_ylim(bottom=1e-4, top=Ymax)
ax.legend(bbox_to_anchor=(1.05, 0.5), loc="center left")
ax.set_xlim(1, 20)

ax.set_xlabel("r [au]")
ax.set_ylabel(r"$\Sigma_d$ normalized to Mdust = 1")

fig.savefig("plot.jpg", bbox_inches="tight", dpi=150)

