# Circular orbit test

## Goal

Check whether the Nbody integrator (Rebound) is properly integrated.

## Rational

We assume rebound to work to machine precision using the IAS15 integrator.
However, we could introduce bugs in our code that destroy this accuracy.
In this test we check that deviations of a circular orbit from the respective sin and cos functions are smaller than the threshold specified in testconfig.yml.