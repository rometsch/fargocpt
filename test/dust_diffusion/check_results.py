import numpy as np
import matplotlib.pyplot as plt
import yaml
from datetime import datetime

from load_dust import get_sigma_dust
from fargocpt import Loader

def test(_):
    datadir = "../../output/tests/dust_diffusion/out"


    fig, axs = plt.subplots(2, 1, dpi=150)
    axs = axs.flatten()
    ax = axs[0]

    toffset = 0

    Ymax = 0
    n = 1
    t = Loader(datadir).snapshots[n]
    sigma_dust, rmid, dr = get_sigma_dust(datadir, n, nbins=101)
    factor = np.sum(sigma_dust*rmid*dr*2*np.pi)
    Y = sigma_dust/factor
    line = ax.plot(rmid, Y, drawstyle="steps-mid", label=f"t = {t:.3f}")

    Ymax = max(Ymax, np.max(Y))


    rref, sigmaref = np.loadtxt("reference_data.txt", unpack=True)

    sigmad_ref = np.interp(rmid, rref, sigmaref)

    dr = np.pad(rmid[1:]-rmid[:-1], (0,1), mode="edge")
    factor = np.sum(sigmad_ref*rmid*dr*2*np.pi)
    reference_data = sigmad_ref/factor



    ax.plot(rmid, reference_data, color="C1", ls="--", label=f"reference")
    Ymax = max(Ymax, np.max(reference_data))
        
    ax.set_yscale("log")
    ax.set_ylim(bottom=1e-6, top=Ymax)
    ax.legend(loc="best")
    ax.set_xlim(1, 20)

    ax.set_xlabel("r [au]")
    ax.set_ylabel(r"$\Sigma_d$ normalized to Mdust = 1")



    diff = np.abs(Y - reference_data)/Ymax

    ax = axs[1]

    ax.plot(rmid, diff)
    ax.set_yscale("log")
    ax.set_ylabel("difference / max value")

    for ax in axs:
        ax.set_xlim(8, 12)


    fig.savefig("plot.jpg", bbox_inches="tight", dpi=150)



    with open("testconfig.yml", "r") as f:
        data = yaml.safe_load(f)
    testname = data["testname"]
    threshold = float(data["threshold"])

    if np.max(diff) < threshold:
        print(f"SUCCESS: {testname}")
    else:
        print(f"FAIL: {testname}, could be statistical, max diff = {np.max(diff)}, theshold = {threshold}")

    with open("test.log", "w") as f:
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"{current_time}", file=f)
        print(f"max diff = {np.max(diff)}, theshold = {threshold}", file=f)


if __name__=="__main__":
    test("dummy")