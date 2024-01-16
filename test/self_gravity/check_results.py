import os
import yaml
import numpy as np
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt
from fargocpt import Loader

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
    filename = f"gr_direct_{Nr}x{Naz}.txt"
    if (os.path.exists(filename)):
        r, gr = np.loadtxt(filename, unpack=True)
        return r, gr
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


    gr_direct = accr[:,1].to_value("cm/s2")
    r_direct = R[:,0].to_value("au")

    # save to file
    np.savetxt(filename, np.array([r_direct, gr_direct]).T, header=f"#Self-gravity acceleration for a {Nr}x{Naz} grid\n#syntax: r_direct gr_direct")

    return r_direct, gr_direct

def test(output_dir):

    ld = Loader(output_dir)

    r_direct, gr_direct = get_gr_direct(ld)    

    # plot difference
    r_code, gr_code = ld.gas.vars2D.avg("a_sg_rad", 0)
    r_code = r_code.to_value("au")
    gr_code = gr_code.to_value("cm/s2")


    fig, axs = plt.subplots(2,1,height_ratios=[3,1], sharex=True, figsize=(6,4))
    fig.subplots_adjust(hspace=0)
    ax = axs[0]
    ax.plot(r_direct, gr_direct, label="direct sum")
    ax.plot(r_code, gr_code, ls="--", label="code")
    ax.set_xlabel("r [au]")
    ax.set_ylabel("$g_r$ [cm/s2]")
    ax.legend()
    ax.axhline(0.0, ls="-", lw=1, color="k", alpha=0.5, zorder=0)

    ax = axs[1]
    diff = np.abs(gr_code/gr_direct - 1)
    ax.plot(r_code, diff, label="relative difference")
    # ax.plot(r_code, np.abs(gr_code-gr_direct)*1e5, label="absolute difference * $10^5$", ls="--", color="C2")
    ax.set_yscale("log")
    ax.set_xlabel("r [au]")
    ax.set_ylabel("relative difference")
    ax.grid(alpha=0.5)
    ax.set_yticks([1e-3, 1e-2, 1e-1, 1])
    # ax.legend(bbox_to_anchor=(1.0, 3.3), loc="upper right")

    kzero = np.argmin(np.abs(gr_code))
    rzero = r_code[kzero]
    for ax in axs:
        ax.axvline(rzero, ls="-", lw=1, color="k", alpha=0.5, zorder=0)
    

    fig.savefig("plot.jpg", dpi=150)
    fig.savefig("plot.pdf", dpi=150)

    threshold = 0.001
    diff_test = diff[r_code > 2]
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
    output_dir = "../../output/tests/self_gravity_solver/out"
    test(output_dir)