import numpy as np
import astropy.units as u
import astropy.constants as const



def vdrift_theo(stokes, r):
    """ Drift speed according to Picogna & Kley 2015 Eq. (C.1) (10.1051/0004-6361/201526921)
    and Nakagawa+1986 Eq. (1.9) (10.1016/0019-1035(86)90121-1). 
    Note that in Eq. (1.9) for eta, there is a '/' missing and r OmegaK^2 needs to be in the denominator.  
    The implementation matches the values in the plot of Picogna & Kley 2015 Fig. C.2 
    (though the mislabeled the unit on the y axis to be cm/s but it must be au/yr)."""
    sigmaslope = -1
    temperatureslope = -1
    Mstar = 1.0*u.solMass
    h = 0.05
    # r = 1*u.au
    vK = np.sqrt(const.G*Mstar/r).decompose()
    eta = h**2 * (sigmaslope + temperatureslope)
    vdrift = eta*vK/(stokes + stokes**-1)
    return vdrift