# MPI simple test

## Goal

Check whether the code runs properly in MPI mode.

## Reason

If MPI barriers are wrongly placed, the code can work with a single process but hang in a multi process scenario.


## Devations close to the star

The code includes corrections for the near field close to the star where the disk scale height is small.
This causes devations from the expected profile in the inner disk.