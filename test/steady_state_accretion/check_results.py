import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import yaml
import os
import sys
import inspect

def test(_):

    # realpath() will make your script run, even if you symlink it :)
    cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    if cmd_folder not in sys.path:
        sys.path.insert(0, cmd_folder)

    # Use this if you want to include modules from a subfolder
    cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(cmd_folder,"../../Tools/read_data")))
    if cmd_subfolder not in sys.path:
        sys.path.insert(0, cmd_subfolder)

    # Info:
    # cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
    # __file__ fails if the script is called in different ways on Windows.
    # __file__ fails if someone does os.chdir() before.
    # sys.argv[0] also fails, because it doesn't not always contains the path.

    # sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'Tools', 'read_data'))
    from read_gasMassFlow1D import gasMassFlow1D
    from read_par_file import read_unit_file
    from read1D import read1D

    out_folder = "../../output/tests/steady_state_accretion/out/"

    units = read_unit_file(out_folder)

    def M_dot_by_visc(radius, surface_density, visc):
        """ Compute the mass flow rate using the viscosity

        Parameters
        ----------
        radius : array
            Radius of the cells
        surface_density : array
            Surface density of the cells
        vis : array
            Viscosity of the cells
        """

        M_dot = 3*np.pi*surface_density * visc
        M_dot = M_dot.decompose().to("solMass/yr")

        return radius, M_dot


    def M_dot_by_vel(radius, surface_density, vr):
        """ Compute the mass flow rate using the radial velocity

        Parameters
        ----------
        radius : array
            Radius of the cells
        surface_density : array
            Surface density of the cells
        vr : array
            Radial velocity of the cells
        """

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



    n = N_max

    _, sigma0 = dens1D_reader.read(0)
    _, vrad0 = vrad1D_reader.read(0)
    vrad0[vrad0 == 0] = np.max(vrad0[vrad0 != 0])
    r_, sigma_ = dens1D_reader.read(n)
    r__, vrad_ = vrad1D_reader.read(n)
    r___, vis_ = vis1D_reader.read(n)


    x_, y_ = M_dot_by_visc(r_, sigma_, vis_)
    x2_, y2_ = M_dot_by_vel(r__, sigma_, vrad_)

    rs, data = mass_flow_reader.read(n)

    data = data.to("solMass/yr")

    Mdot_theo = 1e-8 * u.solMass/u.yr

    diffval_ = np.abs(data[1:-1]) / Mdot_theo - 1
    diff_ = np.abs(y_) / Mdot_theo - 1
    diff2_ = np.abs(y2_) / Mdot_theo - 1

    axs[0].plot(x_, diff_, 's', color='C0' , label='From Viscosity')
    axs[0].plot(x2_, diff2_, '.', color='C2', label='From Velocity')
    axs[0].plot(rs[1:-1], diffval_, label="Simulation " + str(n), color='C1')

    axs[1].plot(r_, (sigma_/sigma0) - 1, '-', label="density " + str(n))
    axs[1].plot(r__, (vrad_/vrad0) - 1, '.', label="vr " + str(n), color='C2')
    axs[1].set_ylabel("relative difference")

    # axs[1].set_ylim(0.01, 100)

    # axs[0].set_ylim(10**-9, 10**-7)
    axs[0].set_ylabel("$\dot{M} / \dot{M}_\mathrm{theo} - 1$")
    axs[0].set_xlabel("au")
    axs[0].legend()
    axs[1].legend()

    for ax in axs:
        ax.set_xlim(20,60)
        ax.set_xscale('log')

    eps = 3e-4
    for ax in axs:
        ax.set_ylim(-eps, eps)

    axs[1].set_xlabel("r [au]")


    ## Check against threshold
    xmin = 20*u.au
    xmax = 60*u.au
    inds = np.logical_and(x_[1:] > xmin, x_[:-1] < xmax)
    maxdiff = np.max(np.abs(diffval_[inds]))

    with open("testconfig.yml", "r") as infile:
        testconfig = yaml.safe_load(infile)

    threshold = testconfig["threshold"]
    testname = testconfig["testname"]
    if (maxdiff < threshold):
        print(f"SUCCESS: {testname}")
    else:
        print(f"FAIL: {testname}")
        
        
    with open("test.log", "w") as outfile:
        print(f"maxdiff = {maxdiff}, threshold = {threshold}", file=outfile)



    fig.tight_layout()
    # plt.savefig('mass_flow.png')
    # plt.show()
    plt.savefig("plot.jpg", dpi=150)

if __name__ == "__main__":
    test("dummy")