# Cold disk test

A power law disk with an ideal equation of state and no heating and cooling is tested for stability.

## Physical setup

The physics setup is as follows:
- ideal equation of state
- no explicit viscosity
- no artificial viscosity, i.e. no shock handling
- no heating from irradiation
- no radiative cooling

## Expected behavior

The expectation is that the disk stays at the initial condition.

## Reason

This test fails when the energy update due to compression heating is performed before the velocity updates from the source terms.