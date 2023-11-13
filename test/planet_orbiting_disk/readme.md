# Planet orbiting disk test

## Goal

Check whether gravity interaction between planets and disk is implemented correctly.

## Rational

For a planet on an orbit it should not matter what the central mass is as long as it is axisymmetric.
Thus, replacing the central star for a disk with the same mass should give the same dynamical result.
This test checks whether this is the case.

## Expectation

The simulation with the first order Euler timestepping should drift rather quickly because the forward Euler integrator is unstable.
The leap frog scheme should give better results.