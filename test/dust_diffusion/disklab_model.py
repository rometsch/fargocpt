import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

try:
    from disklab import DiskRadialModel
    from disklab.natconst import year, au, MS
except ImportError:
    print("Could not import disklab. Make sure it is installed or download it from https://github.com/dullemond/disklab.")
    exit(1)


#-----------------------------------------
# Define physical parameters
#-----------------------------------------

alpha       = 1e-2
tend        = 15_000*year

tstart      = 0.01*year
ntime       = 1_000
nr          = 10_000
h0          = 0.040613
flare       = 0.25

rin         = 0.1*au
rout        = 40*au
r0          = 1*au

rpeak       = 10*au
wpeak       = 0.01*au
grainsize   = 1e-5
sigma0      = 1e2
sigmaslope  = -1
Mstar       = 0.5*MS
Mdisk0      = 1e-4



#-----------------------------------------
# Output handling
#-----------------------------------------

import h5py
import os

def create_outputfile(filename):
    with h5py.File(filename, "w") as dset:
        dset.attrs["label"] = "disklab"
        dset.attrs["unitsystem"] = "cgs"
        dset.create_dataset("time", (nsaves,), dtype=np.float64)
        dset["time"].attrs["unit"] = "s"
        dset.create_dataset("r", (nr,), dtype=np.float64)
        dset["r"].attrs["unit"] = "cm"
        dset.create_dataset("omk", (nr,), dtype=np.float64)
        dset["omk"].attrs["unit"] = "1/s"
        dset.create_dataset("Dd", (nr,), dtype=np.float64)
        dset["Dd"].attrs["unit"] = "cm2/s"
        dset.create_dataset("sigmad", (nsaves,nr), dtype=np.float64)
        dset["sigmad"].attrs["unit"] = "g/cm2"
        dset.create_dataset("n_last_saved", (1,), dtype=np.int64)
        dset["n_last_saved"][0] = -1

def save_constant(r,Dd,omk):
    with h5py.File(filename, "a") as dset:
        dset["r"][:] = r
        dset["Dd"][:] = Dd
        dset["omk"][:] = omk

def save(sigmad,t):
    with h5py.File(filename, "a") as dset:
        n = dset["n_last_saved"][0] + 1
        if n > 0:
            if dset["time"][n-1] == t:
                return
        dset["n_last_saved"][0] = n
        dset["time"][n] = t
        dset["sigmad"][n,:] = sigmad


dirname = "."
filename = f"{dirname}/disklab.hdf5"
nsaves = ntime+1


#-----------------------------------------
# Main loop
#-----------------------------------------


def run(tfinal = tend/year):
    tfinal = tfinal*year

    create_outputfile(filename)

    time   = tstart * (tfinal / tstart)**(np.linspace(0., 1., ntime + 1))
    # time   = np.linspace(tstart, tend, ntime+1)

    d = DiskRadialModel(mstar=Mstar, rin=rin, rout=rout, nr=nr, alpha=alpha)
    d.make_disk_from_powerlaw(sigma0,r0,sigmaslope)

    hr = h0*(d.r/au)**flare # Choice of aspect ratio h_p/r
    d.cs = hr*d.r*d.omk
    d.hp = d.cs/d.omk
    # d.tmid = 2.3*mp*d.cs**2/kk # 2.3 is the mean molecular weight
    d.compute_mass()
    d.compute_rhomid_from_sigma()

    d.add_dust(agrain=grainsize)

    #
    # If requested, compute nu and dmix from alpha
    #
    dust = d.dust[0]
    dust.diskradialmodel.compute_nu()
    if hasattr(dust.diskradialmodel, 'alphamix'):
        dust.dmix = dust.diskradialmodel.alphamix * dust.diskradialmodel.cs * dust.diskradialmodel.cs / dust.omk / dust.diskradialmodel.Sc
    else:
        dust.dmix = dust.diskradialmodel.nu / dust.diskradialmodel.Sc
    dust.dmix[:] *= 1.0 / (1.0 + dust.St**2)


    save_constant(d.r, d.dust[0].dmix, d.omk)


    X = np.exp(-0.5*(d.r - rpeak)**2/ wpeak**2)
    d.dust[0].sigma = X / np.sum(X * d.dsurf) * Mdisk0 * MS

    # d.Sc = 1e10    # Switch off mixing by putting Schmidt number to 'infinity'

    # print("Sc", d.dust[0].diskradialmodel.Sc)
    # d.dust[0].sigma = d.dust[0].get_dust_radial_drift_next_timestep(0.1, fixgas=True)
    # print("dmix", d.dust[0].dmix)
    # exit
    # plotting

    save(d.dust[0].sigma,0)

    for itime in tqdm(range(1, ntime + 1)):
        d.dust[0].sigma = d.dust[0].get_dust_radial_drift_next_timestep(time[itime] - time[itime - 1], fixgas=True)
        save(d.dust[0].sigma,time[itime])

    return d.r/au, d.dust[0].sigma

#-----------------------------------------
# Plotting
#-----------------------------------------

def main():
    from argparse import ArgumentParser

    parse = ArgumentParser()
    parse.add_argument("--plot", action="store_true", help="Plot the last dust profile.")
    opts = parse.parse_args()

    r, sigmad = run()

    if opts.plot:
        plt.plot(r, sigmad)
        plt.xlabel("r [au]")
        plt.ylabel("Sigmad [cgs]")
        plt.show()

if __name__=="__main__":
    main()