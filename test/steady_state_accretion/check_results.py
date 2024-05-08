import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import yaml
from fargocpt import Loader

def test(_):

    out_folder = "../../output/tests/steady_state_accretion/out/"

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
        vr = (vr[1:] + vr[:-1])/2
        M_dot = -2*radius*np.pi*surface_density*vr
        M_dot = M_dot.decompose().to("solMass/yr")

        return radius, M_dot


    ld = Loader(out_folder)
    Nfinal = ld.snapshots[-1]

    fig, axs = plt.subplots(2,1, sharex=True)

    n = Nfinal

    sigma0 = ld.gas.vars1D.avg("Sigma", 0, grid=False)
    vrad0 = ld.gas.vars1D.avg("vrad", 0, grid=False)
    vrad0[vrad0 == 0] = np.max(vrad0[vrad0 != 0])
    r_, sigma_ = ld.gas.vars1D.avg("Sigma", n)
    r__, vrad_ = ld.gas.vars1D.avg("vrad", n)
    r___, vis_ = ld.gas.vars1D.avg("viscosity", n)


    x_, y_ = M_dot_by_visc(r_, sigma_, vis_)
    x2_, y2_ = M_dot_by_vel(r__, sigma_, vrad_)

    rs, data = ld.gas.vars1D.avg("MassFlow", n)

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
    plt.savefig("plot.jpg", dpi=150)

if __name__ == "__main__":
    test("dummy")