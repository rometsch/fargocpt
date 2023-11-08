# 2D test for the FLD module

This test is aimed at the solver of the diffusion equation.
The diffusion equation
$$
\frac{\partial x}{\partial t} = \nabla \cdot \left[ K(x,t) \nabla x \right]
$$
is discretized on the polar grid and the resulting linear equation system is solved using the successive-over-relaxation (SOR) method.

The test consists of an initially guassian profile with a constant diffusion coefficient in the grid to which the analytical solution is known.

The FLD module of the code uses temperature as a variable because it uses the one-temperature approximation in which the radiation energy is neglected against the internal energy density.
In this formulation, the diffusion coefficient depends on the third power of the temperature.
The physics of the FLD module is tested in 1D in another test already.

This test focuses on the diffusion solver in 2D.

To be able to compare to analytical solutions, we take the diffusion coefficient to be constant.

In this test, we solve the equation
$$
\frac{\partial x}{\partial t} = K \Delta x
$$


The FLD module includes extra code for this test.
Only one timestep is performed and the procedure is separated from the usual schedule.
Once the test function is entered, the file `Erad_input.dat` is loaded from within the output directory.
Then the diffusion equation is solved for one timestep and the result is written out to `Erad_output.dat` to the output directory.
This allows to test the diffusion solver independently of the rest of the code.

## Analytical solution

Diffuse the solution to the diffusion problem in 2D at some point in time to a later point in time.
For a two dimensional process this is
$$
x(\vec{r}, t) = \frac{x_0}{4\pi t K} \exp\left( - \frac{|\vec{r} - \vec{r_0}|^2}{4Kt} \right)\,.
$$

The physical variables are taken to be of order one in cgs units.
The parameters are $x_0 = 1\text{K}$, $\vec{r_0} = (1\,\mathrm{cm}, 0)$, $K = 1 \frac{\mathrm{cm}^2}{\mathrm{s}}$.
Note that the unit of $x$ is actually irrelevant for this test, because the value are read just before the diffusion solver is called and written out directly afterwards.

The test is initialized at $t_0 = 10^{-3}\,\mathrm{s}$ and evolved for $\Delta t = 0.9\times10^{-3}\,\mathrm{s}$ until $t = 10^{-2}\,\mathrm{s}$.
This ensures that the spread is still within the grid to avoid boundary effects.

## Numerical setup

The unit system is cgs, with cm, g, seconds, and Kelvin as base units.
The polar grid streches radially from 0.1 cm to 10 cm and is spaced evenly.

The SOR solver is forced to run for 10000 iterations by setting the tolerance to 0.

## Run the test

There is a script `run_auto_test.sh` which sets up the initial conditions and runs the code and checks the deviation of the numerical solution to the analytical solution against a threshold configured in `test_settings.yml`. 
This print a success or fail message suitable for the test suite.

To manually run the test, you can run the following commands.

First run the `run_code.sh` script to setup the initial conditions and run the code to solve the diffusion equation.
Then run `check_differences.py --diff` to display the results and plot the difference.

To check the linear system inversion by the SOR solver, you can call `check_differences.py --diff --solve` which solves the linear system exactly using matrix inversion.
This generally works down to maschine precision give enough iterations for the SOR solver.