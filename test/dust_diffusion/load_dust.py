import numpy as np

from fargocpt import Loader


def get_sigma_dust(outdir, N, nbins=51):

    l = Loader(outdir)
    r = l.particles.get("r", N).value
    counts, interfaces = np.histogram(r, bins=nbins)
    mid = 0.5*(interfaces[1:] + interfaces[:-1])
    dr = interfaces[1:] - interfaces[:-1]
    sigma_dust = counts/(dr*mid*2*np.pi)
    return sigma_dust, mid, dr

