#!/usr/bin/env python3
import argparse

import astropy.constants as const
import astropy.units as u
import numpy as np


def main():
    args = parse_cli_args()

    R0 = u.Quantity(args.R0)
    Sigma0 = u.Quantity(args.Sigma0)
    SigmaSlope = args.SigmaSlope
    T0 = u.Quantity(args.T0)
    Tslope = args.Tslope
    mu = args.mu
    gamma = args.gamma
    Mstar = u.Quantity(args.Mstar)
    outR = u.Quantity(args.outR)

    disk = LocallyIsothermalDisk(
        Mstar=Mstar, gamma=gamma, mean_molweight=mu, R0=outR)

    disk.set_temperature_params(T0, R0, Tslope)
    disk.set_Sigma_params(Sigma0, R0, SigmaSlope)

    disk.temperature_to_aspect_ratio()

    print(disk)


def parse_cli_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("R0", help="Reference radius.")
    parser.add_argument("Sigma0", help="Surface density at reference radius.")
    parser.add_argument("SigmaSlope", type=float,
                        help="Surface density power law exponent.")
    parser.add_argument("T0", help="Temperature at reference radius.")
    parser.add_argument("Tslope", type=float,
                        help="Temperature power law exponent.")
    parser.add_argument("Mstar", help="Stellar mass.")
    parser.add_argument("--mu", default=2.4, type=float,
                        help="Mean molecular weight.")
    parser.add_argument("--gamma", default=1.0, type=float,
                        help="Adiabatic exponent.")
    parser.add_argument("--outR", default="1au",
                        help="Radius at which to print values.")
    args = parser.parse_args()
    return args


class LocallyIsothermalDisk:
    def __init__(self,
                 Mstar=None,
                 gamma=None,
                 mean_molweight=None,
                 alpha=None,
                 R0=1*u.au,
                 name="generic locally isothermal disk"):
        self.R0 = R0
        self.Mstar = Mstar
        self.gamma = gamma
        self.mean_molweight = mean_molweight
        self.alpha = alpha
        self.name = name
        self.Sigma = None
        self.aspect_ratio = None
        self.temperature = None

    def set_Sigma_params(self, Sigma0, R0, slope, exp_cutoff_radius=-1):
        """ set parameters for powerlaw (+ cutoff) parametrization """
        self.Sigma = CutoffPowerLaw(R0,
                                    Sigma0,
                                    slope,
                                    cutoff=exp_cutoff_radius)

    def calc_Sigma(self, R):
        try:
            rv = self.Sigma(R)
        except TypeError:
            raise AttributeError("Can not calculate Sigma for this disk.")
        rv = rv.decompose().to("g/cm2")
        return rv

    def set_aspect_ratio(self, h0, R0, slope):
        """ set parameters for powerlaw parametrization """
        self.aspect_ratio = CutoffPowerLaw(R0, h0, slope)

    def calc_aspect_ratio(self, R):
        try:
            rv = self.aspect_ratio(R)
        except TypeError:
            raise AttributeError(
                "Can not calculate aspect ratio for this disk.")
        rv = rv.decompose().value
        return rv

    def set_temperature_params(self, T0, R0, slope):
        """ set parameters for powerlaw parametrization """
        self.temperature = CutoffPowerLaw(R0, T0, slope)

    def calc_temperature(self, R):
        try:
            rv = self.temperature(R)
        except TypeError:
            raise AttributeError(
                "Can not calculate temperature for this disk.")
        rv = rv.decompose().to("K")
        return rv

    def h_to_T_factor(self):
        Mstar = self.Mstar
        gamma = self.gamma
        mu = self.mean_molweight
        if any([x is None for x in [Mstar, gamma, mu]]):
            raise ValueError("a was not specified")

        factor = mu * 1.008 * u.u / (gamma *
                                     const.k_B) * const.G * Mstar / self.R0
        return factor.decompose()

    def aspect_ratio_to_temperature(self):
        factor = self.h_to_T_factor()
        h0 = self.aspect_ratio(self.R0)
        T0 = factor * h0**2
        T0 = T0.decompose().to("K")
        slope = 2 * self.aspect_ratio.slope - 1
        self.temperature = CutoffPowerLaw(self.R0, T0, slope)

    def temperature_to_aspect_ratio(self):
        factor = self.h_to_T_factor()
        T0 = self.temperature(self.R0)
        h0 = np.sqrt(T0 / factor).decompose().to_value(
            u.dimensionless_unscaled)
        slope = (self.temperature.slope + 1) / 2
        self.aspect_ratio = CutoffPowerLaw(self.R0, h0, slope)

    def calc_mass_accretion_rate(self, R):
        """ calculate the mass accretion rate with Mdot = 3pi Sigma nu and nu = alpha cs H"""
        alpha = self.alpha
        Mstar = self.Mstar
        if alpha is None:
            raise AttributeError("alpha was not specified")
        Sigma = self.Sigma(R)
        hsq = self.aspect_ratio(R)**2
        rv = 3*np.pi*Sigma*alpha*hsq*np.sqrt(const.G*Mstar*R)
        rv = rv.decompose().to("solMass/yr")
        return rv

    def __str__(self):
        rv = "Locally isothermal disk model\n"
        rv += "name = {}\n".format(self.name)
        rv += f"\nValues at R0 = {self.R0}"
        rv += "\nSigma:\n"
        rv += str(self.Sigma(self.R0))
        rv += f"\nslope = {self.Sigma.slope}"
        rv += "\n\ntemperature:\n"
        rv += str(self.temperature(self.R0))
        rv += f"\nslope = {self.temperature.slope}"
        rv += "\n\naspect_ratio:\n"
        rv += str(self.aspect_ratio(self.R0))
        rv += f"\nslope = {self.aspect_ratio.slope}"
        try:
            rv += f"\nmass accretion rate at {self.R0}:\n"
            rv += str(self.calc_mass_accretion_rate(self.R0))
        except AttributeError:
            pass
        return rv


class CutoffPowerLaw:
    def __init__(self, x0, y0, slope, cutoff=None):
        self.x0 = x0
        self.y0 = y0
        self.slope = slope
        self.cutoff = cutoff

    def __call__(self, x):
        rv = self.y0 * (x / self.x0)**self.slope
        if self.cutoff is not None and self.cutoff > 0:
            rv *= np.exp(-x / self.cutoff)
        return rv

    def __str__(self):
        rv = "x0 = {}\n".format(self.x0)
        rv += "y0 = {}\n".format(self.y0)
        rv += "slope = {}\n".format(self.slope)
        return rv


if __name__ == "__main__":
    main()
