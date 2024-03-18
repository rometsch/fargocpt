# Parameters
Following is a list of parameters that can be used to configure physics and numerics of a simulation and the behavior of the program. 

## Active Parameters
 For numerical parameters, the range of choices are abreviated as follows:
- '-' : the range includes negative numbers.
- '0' : the range includes zero.
- '+' : the range includes positive numbers.

| Parameter                | Choices   | Default        | Type   | Unit Support   | Description                                                                                  |
|:-------------------------|:----------|:---------------|:-------|:---------------|:---------------------------------------------------------------------------------------------|
| accretion efficiency     | 0+        | 0              | double | False          | Accretion efficiency factor. Fraction of the mass in the Hill sphere to accrete every orbit. |
| argument of pericenter   | 0-2pi     | 0              | double | False          | Argument of pericenter of the initial orbit in radians.                                      |
| eccentricity             | 0-1       | 0              | double | False          | Eccentricity of the initial orbit.                                                           |
| irradiation ramp-up time | 0+        | 0              | double | True           | Ramp-up time for the irradiation in multiples of the orbital period.                         |
| mass                     | 0+        | none           | double | True           | Inital mass.                                                                                 |
| name                     |           | none           | string | False          | Name of the body.                                                                            |
| radius                   | 0+        | 0.009304813 au | double | True           | Radius of the body (assumed to be a sphere). Is constant in time.                            |
| ramp-up time             | 0+        | 0              | double | False          | Ramp-up time for the mass of the planet. This only affects the gravitational interactions.   |
| semi-major axis          | 0+        | none           | double | True           | Semimajor axis of the initial orbit.                                                         |
| temperature              | 0+        | 0.0 K          | double | True           | Temperature of the body. Is constant in time.                                                |
| trueanomaly              | 0+        | 0              | double | False          | True anomaly of the initial orbit in radians.                                                |
