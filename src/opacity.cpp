#include "opacity.h"
#include "LowTasks.h"
#include "parameters.h"
#include <cassert>
#include <cmath>

namespace opacity
{

double opacity(double density, double temperature)
{
	const double temperatureCGS = temperature * units::temperature.get_code_to_cgs_factor();
	const double densityCGS = density * units::density.get_code_to_cgs_factor();

	double rv;

    switch (parameters::opacity) {
    case parameters::opacity_lin:
	rv = lin(densityCGS, temperatureCGS) * units::opacity.get_cgs_to_code_factor();
	break;

    case parameters::opacity_bell:
	rv = bell(densityCGS, temperatureCGS) * units::opacity.get_cgs_to_code_factor();
	break;

    case parameters::opacity_const_op:
	rv = parameters::kappa_const;
	break;

    case parameters::opacity_simple:
	// this is only used for the temperature test
	rv = parameters::kappa_const * std::pow(temperatureCGS, 2);
	break;

    default:
	die("Invalid opacity!");
	rv = 0;
	break;
    }

	return parameters::kappa_factor * rv;
}

/**
	Opacities after Lin & Papaloizou (1985)
	evaluate kappa (opacity) as functions of density and temperature
	Doug's Opacities for Ice, grains, slightly modified
*/
double lin(double density, double temperature)
{
    // temperature must be positive
    assert(temperature > 0);
    // temperatures above 1.000.000.000 K should not appear
    assert(temperature < 1e9);
    // surface density must be positive
    assert(density > 0);

    const double power1 = 4.44444444e-2;
    const double power2 = 2.381e-2;
    const double power3 = 2.267e-1;

    const double t234 = 1.6e3;
    const double t456 = 5.7e3;
    const double t678 = 2.28e6;

    /* coefficients for opacity laws 1, 2, and 3 in cgs units */
    const double ak1 = 2.e-4;
    const double ak2 = 2.e16;
    const double ak3 = 5.e-3;

    /* coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units */
    const double bk3 = 50.;
    const double bk4 = 2.e-2;
    const double bk5 = 2.e4;
    const double bk6 = 1.e4;
    const double bk7 = 1.5e10;
    const double bk8 = 0.348;

    /* test T against (T_23 * T_34 * T_34)**0.333333333 */
    if (temperature > t234 * std::pow(density, power1)) {
	/* to avoid overflow */
	double ts4 = 1.e-4 * temperature;
	double density13 = std::pow(density, 1.0 / 3.0);
	double density23 = density13 * density13;
	double ts42 = ts4 * ts4;
	double ts44 = ts42 * ts42;
	double ts48 = ts44 * ts44;

	/* test T against (T_45 * T_56)**0.5 */
	if (temperature > t456 * std::pow(density, power2)) {
	    if ((temperature < t678 * std::pow(density, power3)) ||
		(density <= 1e-10)) {
		/* disjoint opacity laws for 5, 6, and 7 */
		double o5 = bk5 * density23 * ts42 * ts4;
		double o6 = bk6 * density13 * ts48 * ts42;
		double o7 = bk7 * density / (ts42 * std::sqrt(ts4));

		/* parameters used for smoothing */
		double o6an = o6 * o6;
		double o7an = o7 * o7;

		/* smoothed and continuous opacity law for regions 5, 6, and 7
		 */
		return std::pow(
		    std::pow(o6an * o7an / (o6an + o7an), 2) +
			std::pow(o5 / (1.0 +
				       std::pow(ts4 / (1.1 * std::pow(density,
								      0.04762)),
						10.0)),
				 4.0),
		    0.25);
	    } else {
		/* disjoint opacity laws for 7 and 8 */
		double o7 = bk7 * density / (ts42 * std::sqrt(ts4));
		double o8 = bk8;

		/* parameters used for smoothing */
		double o7an = o7 * o7;
		double o8an = o8 * o8;

		/* smoothed and continuous opacity law for regions 7 and 8 */
		return std::pow(o7an * o7an + o8an * o8an, 0.25);

		/* no scattering */
		// return bk7*density/(ts42*sqrt(ts4));
	    }
	} else {
	    /*  disjoint opacity laws for 3, 4, and 5 */
	    double o3 = bk3 * ts4;
	    double o4 = bk4 * density23 / (ts48 * ts4);
	    double o5 = bk5 * density23 * ts42 * ts4;
	    /* parameters used for smoothing */
	    double o4an = std::pow(o4, 4);
	    double o3an = std::pow(o3, 4);

	    /* smoothed and continuous opacity law for regions 3, 4, and 5 */
	    return pow((o4an * o3an / (o4an + o3an)) +
			   pow(o5 / (1.0 + 6.561e-5 / ts48), 4),
		       0.25);
	}
    } else {
	/* different powers of temperature */
	double t2 = temperature * temperature;
	double t4 = t2 * t2;
	double t8 = t4 * t4;
	double t10 = t8 * t2;

	/* disjoint opacity laws */
	double o1 = ak1 * t2;
	double o2 = ak2 * temperature / t8;
	double o3 = ak3 * temperature;

	/* parameters used for smoothing */
	double o1an = o1 * o1;
	double o2an = o2 * o2;

	/* smoothed and continuous opacity law for regions 1, 2, and 3 */
	return std::pow(std::pow(o1an * o2an / (o1an + o2an), 2) +
			    std::pow(o3 / (1 + 1.e22 / t10), 4),
			0.25);
    }
}

/**
	Opacities after Bell & Lin (1994)

	OPACITIES UP TO 1-1,000,000 K
	The opacity table is divided into eight regions:
	region 1 is dominated by ice grains
	region 2 is where ice grains melts
	region 3 is dominated by metal grains
	region 4 is where metal grains melts
	region 5 is dominated by molecules
	region 6 is dominated by H-
	region 7 is dominated by Kramer's law
	region 8 is dominated by electron scattering

	The opacity for these eight regions are obtained from fitting
   Bodemheimer's table to power laws of T and rho. A smooth function is
   introduced to bring about the transition between each region

	\param density density
	\param temperature temperature
	\returns opacity
*/
double bell(double density, double temperature)
{
    const double power1 = 2.8369e-2;
    const double power2 = 1.1464e-2;
    const double power3 = 2.2667e-1;

    const double t234 = 1.46e3;
    const double t456 = 4.51e3;
    const double t678 = 2.37e6;

    /* coefficients for opacity laws 1, 2, and 3 in cgs units */
    const double ak1 = 2.e-4;
    const double ak2 = 2.e16;
    const double ak3 = 0.1e0;

    /* coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units */
    const double bk3 = 10.;
    const double bk4 = 2.e-15;
    const double bk5 = 1e4;
    const double bk6 = 1e4;
    const double bk7 = 1.5e10;
    const double bk8 = 0.348;

    if (temperature < 1.0) {
	temperature = 10.0;
    }

    if (temperature > t234 * std::pow(density, power1)) {
	/* to avoid overflow */
	double ts4 = 1.e-4 * temperature;
	double density13 = std::pow(density, 1.0 / 3.0);
	double density23 = density13 * density13;
	double ts42 = ts4 * ts4;
	double ts44 = ts42 * ts42;
	double ts48 = ts44 * ts44;

	/* test T against (T_45 * T_56)**0.5 */
	if (temperature > t456 * std::pow(density, power2)) {
	    /* test T against (T67 * T78)**.5 */
	    if ((temperature < t678 * std::pow(density, power3)) ||
		(((density <= 1e10) && (temperature < 1e4)))) {
		/* disjoint opacity laws for 5, 6, and 7 */
		double o5 = bk5 * density23 * ts42 * ts4;
		double o6 = bk6 * density13 * ts48 * ts42;
		double o7 = bk7 * density / (ts42 * std::sqrt(ts4));

		/* parameters used for smoothing */
		double o6an = o6 * o6;
		double o7an = o7 * o7;

		/* smoothed and continuous opacity law for regions 5, 6, and 7
		 */
		return std::pow(
		    std::pow((o6an * o7an / (o6an + o7an)), 2) +
			std::pow(
			    (o5 /
			     (1 + std::pow((ts4 /
					    (1.1 * std::pow(density, 0.04762))),
					   10))),
			    4),
		    0.25);
	    } else {
		/* disjoint opacity laws for 7 and 8 */
		double o7 = bk7 * density / (ts42 * std::sqrt(ts4));
		double o8 = bk8;

		/* parameters used for smoothing */
		double o7an = o7 * o7;
		double o8an = o8 * o8;

		/* smoothed and continuous opacity law for regions 7 and 8 */
		return std::pow(o7an * o7an + o8an * o8an, 0.25);
	    }
	} else {
	    /* disjoint opacity laws for 3, 4, and 5 */
	    double o3 = bk3 * std::sqrt(ts4);
	    double o4 = bk4 * density / (ts48 * ts48 * ts48);
	    double o5 = bk5 * density23 * ts42 * ts4;

	    /* parameters used for smoothing */
	    double o4an = std::pow(o4, 4);
	    double o3an = std::pow(o3, 4);

	    /* smoothed and continuous opacity law for regions 3, 4, and 5 */
	    return std::pow(
		(o4an * o3an / (o4an + o3an)) +
		    std::pow(o5 / (1 + 6.561e-5 / ts48 * 1e2 * density23), 4),
		0.25);
	}
    } else {
	/* different powers of temperature */
	double t2 = temperature * temperature;
	double t4 = t2 * t2;
	double t8 = t4 * t4;
	double t10 = t8 * t2;

	/* disjoint opacity laws */
	double o1 = ak1 * t2;
	double o2 = ak2 * temperature / t8;
	double o3 = ak3 * std::sqrt(temperature);

	/* parameters used for smoothing */
	double o1an = o1 * o1;
	double o2an = o2 * o2;

	/* smoothed and continuous opacity law for regions 1, 2, and 3 */
	return std::pow(std::pow(o1an * o2an / (o1an + o2an), 2) +
			    std::pow(o3 / (1 + 1.e22 / t10), 4),
			0.25);
    }
}

} // namespace opacity
