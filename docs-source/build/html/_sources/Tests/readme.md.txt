# Dust Diffusion test

Charnoz et al. 2011 (DOI:10.1088/0004-637X/737/1/33) perform a dust diffusion test in a 2D disk (Sect. 3.3 Radial Diffusion in the Midplane) to test
diffusion in the midplane.

This is a reproduction of their test.

## particle properties

particle size 0.1 micro m resulting in Stokes number of 1e-7 to 1e-8 between 1 and 10 au.

10 000 particles released at 10 au.

## Gas disk model

$\rho_g(r,z) = \rho_0 r^p \exp(-\frac{z^2}{2H(r)})$
with
$H(r) = H_0 r^\frac{q+3}{2}$
, thus $T \propto r^q$.

Surface density profile has an exponent of $s = p + \frac{q+3}{2}$

Here, $p = -2.25$ and $q = -0.5$, thus $s = -1$.

The flaring parameter is $\frac{-0.5 + 1}{2} = \frac{0.5}{2} = 0.25$

At 1 au, $\Sigma = 200 \frac{kg}{m^2} = 20 \frac{g}{cm^2}$ and $T=204 K$, resulting in $h = 0.040613$ with the values below.
At 5.2 au, $\Sigma = 3.846 \frac{g}{cm^2}$ and $h = 0.061328$.

Central star mass is 0.5 solar masses.

The model is based upon [Brauer et al. 2007](aanda.org/articles/aa/abs/2008/12/aa7759-07/aa7759-07.html)
who use a mean molecular weight of 2.3. They use the internal isothermal sound speed $c_s = \sqrt{\frac{k_B T}{\mu}}$ for $H = \frac{c_s}{\Omega_K}$, thus the adiabatic exponent is $\gamma=1$.

Gas not moving radially, $v_r = 0$
viscosity modelled with $\alpha=0.01$ which is then the only source of radial particle movement.



System evolved for 15 000 yr with snapshots at 2500, 5000, 10 000, and 15 000 yr.

For the Fargo system, turn off disk evolution to keep the background fixed.

### Random number generator for stochastic kicks

Sampling from the standard normal distribution is done using the [Ziggurat algorithm](https://en.wikipedia.org/wiki/Ziggurat_algorithm) using a [C++ implementation](https://github.com/snsinfu/cxx-ziggurat).

The random number generation is done using the [JSF algorithm](http://burtleburtle.net/bob/rand/smallprng.html) using [this implementation](https://gist.github.com/imneme/85cff47d4bad8de6bdeb671f9c76c814).

## Results

### Initial test
![Distribution of the dust particles at different times.](particle_dist_diffusion_test_20-11-16.svg)

To be compared to Charnoz et al. 2011.
Here, there is an error in the formula for the variance of the distibution.
There is a square root missing making the variance way too low resulting in too narrow distributions.

The net inward drift due to dust diffusion is clearly visible.
Its magnitude of 5 au per 10 kyr is quite large! 

### Corrected variance calculation and adaptive time stepping

![Distribution of the dust particles at different times with corrected variance](particle_dist_diffusion_test_adaptive_newrand_20-11-20.svg)

There are two kinks appearing.
This could be due to the vanishing time correlation of the random numbers [Laibe et al. 2020](https://academic.oup.com/mnras/article/494/4/5134/5829884).
But it appears at the CPU domain boundaries, so its probably a communication error.

Two tests for the communication error hypothesis.

1. Run on a single core to see whether the issue remains: out-adaptive-newrand-single.

![Particle location histogram on single core](particle-distribution-singlecore.svg)

2. Call the move function after performing dust diffusion, at the very end of integrate: out-adaptive-newrand-movelast.

![Multicore test with move called after diffusion](particle-distribution-multicore-move-after-diff.svg)

The particle location histogram is qualitatively similar to the Charnoz et al. 2011.
Thus, physics is modelled correctly.

However, something goes wrong in the multicore run at the domain boundaries.

The fields for the gas diffusion coefficient and the radial derivative of the density don't differ across runs with different numbers of cores.
Thus, the problem must be related to the movement of particles across domains.

Removing the resizing of the particles vector in the move functions solves the issue of particles accumulating
at the boundaries.
I do not know why this happened at the moment.

Now, the runs with 1, 2, and 3 cores all yield the same result.

![Comparison of the dust diffusions checks for 1-3 cores.](particle-distribution-core-comparison.svg)