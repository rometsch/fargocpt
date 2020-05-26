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

out_folder = "./out_damp_surface_only/"
out_folder = "./out/"
setup  = read_par_file("./fargo.par")
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


mass_flow_reader = read1D(out_folder,"gasMassFlow")
dens1D_reader = read1D(out_folder, "gasdens")
vrad1D_reader = read1D(out_folder, "gasvrad")
vis1D_reader = read1D(out_folder, "gasviscosity")

N_output = np.loadtxt(out_folder + "timeCoarse.dat", usecols=(0,1), dtype=np.int)
radii = np.loadtxt(out_folder + "used_rad.dat")
N_DT = N_output[:,1]
N_output = np.array([0, 30, 60, 90, 140])
N_output = 10+np.array([0, 3, 6, 9])

fig, ax = plt.subplots()

ax.set_yscale('log')
ax.set_xscale('log')

quantities = np.loadtxt("out/Quantities.dat", skiprows=24)

mass = quantities[:,3]



colors = ['blue', 'red', 'gold', 'green', 'black']
ind = 0
for n in N_output:

    print(mass[N_DT[n]])
    r_, sigma_ = dens1D_reader.read(n)
    r__, vrad_ = vrad1D_reader.read(n)
    r___, vis_ = vis1D_reader.read(n)

    x_, y_ = M_dot(r_, sigma_, vis_)
    x2_, y2_ = M_dot2(r__, sigma_, vrad_)

    rs, data = mass_flow_reader.read(n)
    data = data.to("solMass/yr")
    ax.plot(rs[1:-1], np.abs(data[1:-1]), color=colors[ind], label="Simulation " + str(n))
    ax.plot(x_, np.abs(y_), 's', color=colors[ind], label="Theorie " + str(n))
    #ax.plot(r_, np.abs(sigma_), 's',color=colors[ind], label="density" + str(n))
    #ax.plot(r__, np.abs(vrad_), '.',color=colors[ind], label="vr" + str(n))
    ax.plot(x2_, np.abs(y2_), '.',color=colors[ind], label="Theorie2 " + str(n))
    ind += 1



ax.set_ylabel("M_sol/yr")
ax.set_xlabel("au")
plt.legend()


plt.savefig('mass_flow.png')
plt.show()

