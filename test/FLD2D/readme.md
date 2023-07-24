# 2D test for the FLD module

This test is aimed at the solver of the diffusion equation.
The diffusion equation
$$
\frac{\partial x}{\partial t} = \nabla \cdot \left[ K(x,t) \nabla x \right]
$$
is discretized on the polar grid and the resulting linear equation system is solved using the successive-over-relaxation (SOR) method.

The test consists of a diffusion of point source with a constant diffusion coefficient in the grid to which the analytical solution is known.

The FLD module of the code uses temperature as a variable because it uses the one-temperature approximation in which the radiation energy is neglected against the internal energy density.
In this formulation, the diffusion coefficient depends on the third power of the temperature.

To test the solution of the diffusion equation, we instead solve the FLD problem based on the radiation energy alone with a constant diffusion coefficient.
$$
\frac{\partial E_\mathrm{rad}}{\partial t} = \nabla \cdot \left[ \frac{\lambda c}{\rho \kappa} \nabla E_\mathrm{rad}\right]
$$
Thus, temperature takes the role of radiation energy density.

The FLD module includes extra code for this test.
Only one timestep is performed and the proceedure is separated from the usual schedule.
Once the test function is entered, the file `Erad_input.dat` is loaded from within the output directory.
Then the diffusion equation is solved for one timestep and the result is written out to `Erad_output.dat` to the output directory.
This allows to test the diffusion solver independently of the rest of the code.

## Physical problem description

Diffuse a delta function with a constant diffusion coefficient.
The analytical solution to this process is known.
For a two dimensional process this is
$$
E_\mathrm{rad}(\vec{r}, t) = \frac{E_0}{4\pi t K} \exp\left( - \frac{|\vec{r} - \vec{r_0}|^2}{4Kt} \right)\,.
$$

The parameters are $E_0 = 100\,\mathrm{erg}$, $\vec{r_0} = (5\,\mathrm{cm}, 0)$, $K=\frac{\lambda c}{\rho \kappa} \approx 10^{10} \frac{\mathrm{cm}^2}{\mathrm{s}}$ with $\lambda = 1/3$, $\rho = 1 \mathrm{g}/\mathrm{cm}^3$ and $\kappa = 1 \mathrm{cm}^2/\mathrm{s}$.

The test is initialized at $t_0 = 10^{-11}\,\mathrm{s}$ and evolved for $\Delta t = 10^{-10}\,\mathrm{s}$.
This ensures that the spread is still within the grid to avoid boundary effects.

## Numerical setup

The unit system is cgs, with cm, g, seconds, and Kelvin as base units.
The polar grid streches radially from 0.1 cm to 10 cm and is spaced evenly.

The SOR solver is forced to run for 5000 iterations by setting the tolerance to 0.

## Run the test

First run the `run_code.sh` script to setup the initial conditions and run the code to solve the diffusion equation.
Then run `check_differences.py --diff` to display the results and plot the difference.

To check the linear system inversion by the SOR solver, you can call `check_differences.py --diff --solve` which solves the linear system exactly using matrix inversion.
This generally works down to maschine precision give enough iterations for the SOR solver.