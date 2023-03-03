import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u

import os
import sys
import inspect

# realpath() will make your script run, even if you symlink it :)
cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

# Use this if you want to include modules from a subfolder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../../read_data")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

# Info:
# cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
# __file__ fails if the script is called in different ways on Windows.
# __file__ fails if someone does os.chdir() before.
# sys.argv[0] also fails, because it doesn't not always contains the path.

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'Tools', 'read_data'))
from read_gasMassFlow1D import gasMassFlow1D
from read_par_file import read_par_file, read_unit_file
from read1D import read1D

out_folder = "./out/"

units = read_unit_file(out_folder)
def M_dot(radius, surface_density, vis):

    M_dot = 3*np.pi*surface_density * vis
    M_dot = M_dot.decompose().to("solMass/yr")

    return radius, M_dot

def M_dot2(radius, surface_density, vr):

    radius = (radius[1:] + radius[:-1])/2
    M_dot = -2*radius*np.pi*surface_density*(vr[1:]+ vr[:-1])/2
    M_dot = M_dot.decompose().to("solMass/yr")

    return radius, M_dot


mass_flow_reader = read1D(out_folder,"MassFlow")
dens1D_reader = read1D(out_folder, "Sigma")
vrad1D_reader = read1D(out_folder, "vrad")
vis1D_reader = read1D(out_folder, "viscosity")

N_output = np.loadtxt(out_folder + "/monitor/timeMonitor.dat", usecols=(0,1), dtype=int)
radii = np.loadtxt(out_folder + "used_rad.dat")
N_DT = N_output[:,1]
N_max = N_output[-1][0]
N_output = np.array([0, int(N_max*0.333), int(N_max*0.667), N_max])

fig, axs = plt.subplots(2,1, sharex=True)

for ax in axs[:-0]:
    ax.set_yscale('log')
    ax.set_xscale('log')



ind = 0
for n in N_output:

    _, sigma0 = dens1D_reader.read(0)
    _, vrad0 = vrad1D_reader.read(0)
    vrad0[vrad0 == 0] = np.max(vrad0[vrad0 != 0])
    r_, sigma_ = dens1D_reader.read(n)
    r__, vrad_ = vrad1D_reader.read(n)
    r___, vis_ = vis1D_reader.read(n)

    x_, y_ = M_dot(r_, sigma_, vis_)
    x2_, y2_ = M_dot2(r__, sigma_, vrad_)

    rs, data = mass_flow_reader.read(n)
    if n == 0:
        data[data == 0] = np.mean(y_)
    data = data.to("solMass/yr")
    axs[0].plot(rs[1:-1], np.abs(data[1:-1]), label="Simulation " + str(n), color='C' + str(ind))
    axs[0].plot(x_, np.abs(y_), 's', color='C' + str(ind), label='From Viscosity')
    #ax.plot(r_, np.abs(sigma_), 's',label="density" + str(n))
    #ax.plot(r__, np.abs(vrad_), '.',label="vr" + str(n))
    axs[0].plot(x2_, np.abs(y2_), '.', color='C' + str(ind), label='From Velocity')

    axs[1].plot(r_, (sigma_/sigma0), '-', label="density" + str(n))
    axs[1].plot(r__, (vrad_/vrad0), '.', label="vr" + str(n), color='C' + str(ind))

    ind += 1

# axs[1].set_ylim(0.01, 100)

# axs[0].set_ylim(10**-9, 10**-7)
axs[0].set_ylabel("M_sol/yr")
axs[0].set_xlabel("au")
axs[0].legend()
axs[1].legend()


# plt.savefig('mass_flow.png')
plt.show()

