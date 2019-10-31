#include "opacity.h"
#include "LowTasks.h"
#include "parameters.h"
#include "util.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

namespace opacity
{

double opacity(double density, double temperature)
{
    switch (parameters::opacity) {
    case parameters::opacity_lin:
	return lin(density, temperature);
	break;

    case parameters::opacity_bell:
	return bell(density, temperature);
	break;

    case parameters::opacity_zhu:
	return zhu(density, temperature);
	break;

    case parameters::opacity_kramers:
	return kramers(density, temperature);
	break;

    case parameters::opacity_const_op:
	return parameters::kappa_const;
	break;

    default:
	die("Invalid opacity!");
	return 0;
	break;
    }
}

/**
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
    if (temperature > t234 * pow(density, power1)) {
	/* to avoid overflow */
	double ts4 = 1.e-4 * temperature;
	double density13 = pow(density, 1.0 / 3.0);
	double density23 = density13 * density13;
	double ts42 = ts4 * ts4;
	double ts44 = ts42 * ts42;
	double ts48 = ts44 * ts44;

	/* test T against (T_45 * T_56)**0.5 */
	if (temperature > t456 * pow(density, power2)) {
	    if ((temperature < t678 * pow(density, power3)) ||
		(density <= 1e-10)) {
		/* disjoint opacity laws for 5, 6, and 7 */
		double o5 = bk5 * density23 * ts42 * ts4;
		double o6 = bk6 * density13 * ts48 * ts42;
		double o7 = bk7 * density / (ts42 * sqrt(ts4));

		/* parameters used for smoothing */
		double o6an = o6 * o6;
		double o7an = o7 * o7;

		/* smoothed and continuous opacity law for
		 * regions 5, 6, and 7 */
		return pow(
		    pow(o6an * o7an / (o6an + o7an), 2) +
			pow(o5 / (1.0 + pow(ts4 / (1.1 * pow(density, 0.04762)),
					    10.0)),
			    4.0),
		    0.25);
	    } else {
		/* disjoint opacity laws for 7 and 8 */
		double o7 = bk7 * density / (ts42 * sqrt(ts4));
		double o8 = bk8;

		/* parameters used for smoothing */
		double o7an = o7 * o7;
		double o8an = o8 * o8;

		/* smoothed and continuous opacity law for
		 * regions 7 and 8 */
		return pow(o7an * o7an + o8an * o8an, 0.25);

		/* no scattering */
		// return bk7*density/(ts42*sqrt(ts4));
	    }
	} else {
	    /*  disjoint opacity laws for 3, 4, and 5 */
	    double o3 = bk3 * ts4;
	    double o4 = bk4 * density23 / (ts48 * ts4);
	    double o5 = bk5 * density23 * ts42 * ts4;
	    /* parameters used for smoothing */
	    double o4an = pow4(o4);
	    double o3an = pow4(o3);

	    /* smoothed and continuous opacity law for regions 3, 4,
	     * and 5 */
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

	/* smoothed and continuous opacity law for regions 1, 2, and 3
	 */
	return pow(pow(o1an * o2an / (o1an + o2an), 2) +
		       pow(o3 / (1 + 1.e22 / t10), 4),
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

    if (temperature < 1)
	temperature = 10.0;

    if (temperature > t234 * pow(density, power1)) {
	/* to avoid overflow */
	double ts4 = 1.e-4 * temperature;
	double density13 = pow(density, 1.0 / 3.0);
	double density23 = density13 * density13;
	double ts42 = ts4 * ts4;
	double ts44 = ts42 * ts42;
	double ts48 = ts44 * ts44;

	/* test T against (T_45 * T_56)**0.5 */
	if (temperature > t456 * pow(density, power2)) {
	    /* test T against (T67 * T78)**.5 */
	    if ((temperature < t678 * pow(density, power3)) ||
		(((density <= 1e10) && (temperature < 1e4)))) {
		/* disjoint opacity laws for 5, 6, and 7 */
		double o5 = bk5 * density23 * ts42 * ts4;
		double o6 = bk6 * density13 * ts48 * ts42;
		double o7 = bk7 * density / (ts42 * sqrt(ts4));

		/* parameters used for smoothing */
		double o6an = o6 * o6;
		double o7an = o7 * o7;

		/* smoothed and continuous opacity law for
		 * regions 5, 6, and 7 */
		return pow(
		    pow((o6an * o7an / (o6an + o7an)), 2) +
			pow((o5 /
			     (1 +
			      pow((ts4 / (1.1 * pow(density, 0.04762))), 10))),
			    4),
		    0.25);
	    } else {
		/* disjoint opacity laws for 7 and 8 */
		double o7 = bk7 * density / (ts42 * sqrt(ts4));
		double o8 = bk8;

		/* parameters used for smoothing */
		double o7an = o7 * o7;
		double o8an = o8 * o8;

		/* smoothed and continuous opacity law for
		 * regions 7 and 8 */
		return pow(o7an * o7an + o8an * o8an, 0.25);
	    }
	} else {
	    /* disjoint opacity laws for 3, 4, and 5 */
	    double o3 = bk3 * sqrt(ts4);
	    double o4 = bk4 * density / (ts48 * ts48 * ts48);
	    double o5 = bk5 * density23 * ts42 * ts4;

	    /* parameters used for smoothing */
	    double o4an = pow(o4, 4);
	    double o3an = pow(o3, 4);

	    /* smoothed and continuous opacity law for regions 3, 4,
	     * and 5 */
	    return pow((o4an * o3an / (o4an + o3an)) +
			   pow(o5 / (1 + 6.561e-5 / ts48 * 1e2 * density23), 4),
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
	double o3 = ak3 * sqrt(temperature);

	/* parameters used for smoothing */
	double o1an = o1 * o1;
	double o2an = o2 * o2;

	/* smoothed and continuous opacity law for regions 1, 2, and 3
	 */
	return pow(pow(o1an * o2an / (o1an + o2an), 2) +
		       pow(o3 / (1 + 1.e22 / t10), 4),
		   0.25);
    }
}

/**
	Opacities after Zhu et al. (2012)
	based on http://www.astro.princeton.edu/~zhzhu/opaczhuice.c
*/
double zhu(double density, double temperature)
{
    double xlop, xlp, xlt, rosstabd, pre;

    pre = density * temperature * 8.314472 * 1.e7 / 2.4;
    if (pre < 0. || temperature < 0.) {
	fprintf(stderr, "error: pre or T negative\n");
	fprintf(stderr, "pre: %g T: %g\n", pre, temperature);
    }
    if (pre == 0 || temperature == 0)
	return (1.);
    xlp = log10(pre);
    xlt = log10(temperature);

    if (xlt < 2.23567 + 0.01899 * (xlp - 5.)) {
	xlop = 1.5 * (xlt - 1.16331) - 0.736364;
    } else if (xlt < 2.30713 + 0.01899 * (xlp - 5.)) {
	xlop = -3.53154212 * xlt + 8.767726 -
	       (7.24786 - 8.767726) * (xlp - 5.) / 16.;
    } else if (xlt < 2.79055) {
	xlop = 1.5 * (xlt - 2.30713) + 0.62;
    } else if (xlt < 2.96931) {
	xlop = -5.832 * xlt + 17.7;
    } else if (xlt < 3.29105 + (3.29105 - 3.07651) * (xlp - 5.) / 8.) {
	xlop = 2.129 * xlt - 5.9398;
    } else if (xlt < 3.08 + 0.028084 * (xlp + 4)) {
	xlop = 129.88071 - 42.98075 * xlt +
	       (142.996475 - 129.88071) * 0.1 * (xlp + 4);
    } else if (xlt < 3.28 + xlp / 4. * 0.12) {
	xlop = -15.0125 + 4.0625 * xlt;
    } else if (xlt < 3.41 + 0.03328 * xlp / 4.) {
	xlop = 58.9294 - 18.4808 * xlt + (61.6346 - 58.9294) * xlp / 4.;
    } else if (xlt < 3.76 + (xlp - 4) / 2. * 0.03) {
	xlop = -12.002 + 2.90477 * xlt + (xlp - 4) / 4. * (13.9953 - 12.002);
    } else if (xlt < 4.07 + (xlp - 4) / 2. * 0.08) {
	xlop = -39.4077 + 10.1935 * xlt + (xlp - 4) / 2. * (40.1719 - 39.4077);
    } else if (xlt < 5.3715 + (xlp - 6) / 2. * 0.5594) {
	xlop = 17.5935 - 3.3647 * xlt + (xlp - 6) / 2. * (17.5935 - 15.7376);
    } else {
	xlop = -0.48;
    }

    if (xlop < 3.586 * xlt - 16.85 && xlt < 4.) {
	xlop = 3.586 * xlt - 16.85;
    }

    rosstabd = pow(10., xlop);

    return rosstabd;
}

double kramers(double density, double temperature)
{
    double kappa_kramer, thomson_scattering;

    thomson_scattering = 0.335;
    kappa_kramer = 5.e24 * density * pow(temperature, -3.5);

    return thomson_scattering + kappa_kramer;
}

}; // namespace opacity
