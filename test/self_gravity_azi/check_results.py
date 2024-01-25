import os
import yaml
import numpy as np
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt
from fargocpt import Loader
import pickle

def smoothing_length_sq(r, rp, h):
    """ Squared smoothing length for the SG interaction following Moldenhauer 2018.
    
    The factors in the polinomials are for a grid with rout/rin = 12.5.

    Parameters
    ----------
    r : float
        Radial position cell 1
    rp : float
        Radial position cell 2
    h : float
        Aspect ratio of the disk.

    Returns
    -------
    float
        Squared smoothing length.    
    """
    chi = 0.6472*h - 0.7543*h**2
    lam = 0.4571*h + 0.6737*np.sqrt(h)
    eps_sq = lam**2 * (r - rp)**2 + chi**2 * r*rp
    return eps_sq

def get_gr_direct(ld: Loader):
    """ Return the cached values or calculate the radial acceleration due to SG by direct summation."""
    Nr = ld.gas.vars2D.grid.Nrad
    Naz = ld.gas.vars2D.grid.Naz
    filename = f"gr_direct_{Nr}x{Naz}"
    if (os.path.exists(filename+".pkl")):
        with open(filename+".pkl", "rb") as f:
            r, gr, gaz = pickle.load(f)
        return r, gr, gaz
    # compute the values

    # get values from output
    R, PHI, Sigma = ld.gas.vars2D.get("Sigma", 0, grid_for_plot=False)

    X = R*np.cos(PHI)
    Y = R*np.sin(PHI)
    A = ld.gas.grid.Agrid

    # Calculate SG acceleration in cylindrical coordinates
    accx = np.zeros(Sigma.shape)*u.cm/u.s**2
    accy = np.zeros(Sigma.shape)*u.cm/u.s**2

    h0 = ld.params["AspectRatio"]
    G = const.G

    Nr = Sigma.shape[0]
    Naz = Sigma.shape[1]

    # loop over all cells to compute acceleration there
    for n in range(Nr):
        for k in range(Naz):
            # now use numpy arrays for the calculation

            r1 = R[n,k]

            eps_sq = smoothing_length_sq(r1, R, h0)

            x1 = X[n,k]
            y1 = Y[n,k]

            dx = x1 - X
            dy = y1 - Y
            dsq = dx**2 + dy**2
            d = np.sqrt(dsq)
            
            accx_full = - G * A * Sigma * dx / (d**2 + eps_sq)**1.5
            accy_full = - G * A * Sigma * dy / (d**2 + eps_sq)**1.5

            accx[n,k] = np.sum(accx_full)
            accy[n,k] = np.sum(accy_full)

    # convert to radial acceleration
    accr = accx * np.cos(PHI) + accy * np.sin(PHI)
    accaz = -accx * np.sin(PHI) + accy * np.cos(PHI)


    gr_direct = accr.to_value("cm/s2")
    gaz_direct = accaz.to_value("cm/s2")
    r_direct = R.to_value("au")

    # save to file
    # np.savetxt(filename, np.array([r_direct, gr_direct, gaz_direct]).T, header=f"#Self-gravity acceleration for a {Nr}x{Naz} grid\n#syntax: r_direct gr_direct gaz_direct")
    with open(filename+".pkl", "wb") as f:
        pickle.dump(np.array([r_direct, gr_direct, gaz_direct]), f)

    return r_direct, gr_direct, gaz_direct

def test(output_dir, format="jpg"):

    ld = Loader(output_dir)

    r_direct, gr_direct, gaz_direct = get_gr_direct(ld)

    # plot difference
    R_code, PHI_code, gr_code = ld.gas.vars2D.get("a_sg_rad", 0, grid_for_plot=True)
    R_code, PHI_code, gaz_code = ld.gas.vars2D.get("a_sg_azi", 0, grid_for_plot=True)

    
    r_code = R_code[:,0]
    r_code = r_code.to_value("au")
    gr_code = gr_code.to_value("cm/s2")
    gaz_code = gaz_code.to_value("cm/s2")

    diff_r = np.abs(gr_code - gr_direct)
    rel_diff_r = np.abs(gr_code/gr_direct - 1)
    
    diff_az = np.abs(gaz_direct - gaz_code)
    rel_diff_az = np.abs(gaz_code/gaz_direct - 1)
    
    
    # print(gaz_code - gaz_direct)
    print(gaz_code)

    X_code = R_code*np.cos(PHI_code)
    Y_code = R_code*np.sin(PHI_code)
    
    X_code = X_code.to_value("au")
    Y_code = Y_code.to_value("au")
    
    
    
    Sigma = ld.gas.vars2D.get("Sigma", 0, grid=False)
    # print(Sigma)
    fig, ax = plt.subplots(1,1)

    im = ax.pcolormesh(X_code, Y_code, Sigma.to_value("g/cm2"), cmap="viridis")
    cbar = fig.colorbar(im)
    ax.set_aspect("equal")
    ax.set_xlabel("X [au]")
    ax.set_ylabel("Y [au]")
    cbar.set_label(r"$\Sigma$ [g/cm$^2$]")
    
    fig.savefig(f"Sigma.{format}", dpi=150)
    
    
    fig, ax = plt.subplots(1,1)
    im = ax.pcolormesh(X_code, Y_code, gaz_direct, cmap="bwr")
    cbar = fig.colorbar(im)
    ax.set_aspect("equal")
    ax.set_xlabel("X [au]")
    ax.set_ylabel("Y [au]")
    cbar.set_label(r"$g_\phi$ [cm/s$^2$]")
    fig.savefig(f"gaz.{format}", dpi=150)
    
    rs = ld.gas.grid.radc.to_value("au")
    phis = ld.gas.grid.phic
    
    krad = np.argmin(np.abs(r_code - 4))

    fig, axs = plt.subplots(2,1,height_ratios=[3,1], sharex=True, figsize=(6,4))
    fig.subplots_adjust(hspace=0)
    ax = axs[0]
    ax.plot(phis, gaz_direct[krad,:], label="direct sum", lw=2)
    ax.plot(phis, gaz_code[krad,:], ls="--", label="code", lw=2)
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel("$g_\phi$ [cm/s2]")
    ax.legend()
    ax.axhline(0.0, ls="-", lw=1, color="k", alpha=0.5, zorder=0)
    ax.set_xlim(0, 2*np.pi)

    ax = axs[1]
    ax.plot(phis, rel_diff_az[krad,:], label="rel. diff. azi", lw=2)#, ls="none", marker=".")
    ax.plot(phis, diff_az[krad,:]*1e6, label=r"abs. diff. azi $\times 10^6$", ls="-", lw=2)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel("(relative) difference")
    # ax.grid(alpha=0.5)
    ax.set_yticks([1e-6, 1e-4, 1e-2])
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.01))
    ax.set_xlim(0, 2*np.pi)

    # kzero = np.argmin(np.abs(gr_code))
    # rzero = r_code[kzero]
    # for ax in axs:
    #     ax.axvline(rzero, ls="-", lw=1, color="k", alpha=0.5, zorder=0)
    

    fig.savefig(f"plot.{format}", dpi=300, bbox_inches="tight")

    threshold = 0.0012
    diff_test = diff_r[rs > 2]
    max_diff = np.max(diff_test)

    with open("testconfig.yml", "r") as f:
        testconfig = yaml.safe_load(f)

    testname = testconfig["testname"]
    threshold = testconfig["threshold"]

    pass_test = max_diff < threshold

    with open("test.log", "w") as f:
        print(f"Test name: {testname}", file=f)
        print(f"Max rel diff outside 2 au = {max_diff:.2e}", file=f)
        print(f"Threshold: {threshold}", file=f)
        print(f"Pass test: {pass_test}", file=f)

    if pass_test:
        print(f"SUCCESS: {testname}")
    else:
        print(f"FAIL: {testname}")


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", default="jpg")
    opts = parser.parse_args()
    
    output_dir = "../../output/tests/self_gravity_solver_azi/out"
    test(output_dir, format=opts.format)