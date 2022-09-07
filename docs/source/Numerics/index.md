# Physics and numerics

This part of the documentation describes the numerics choices and algorithms used in the FargoCPT code.

## Topics to include

- hydro simulations
  - transport step with reference to FARGO paper
  - grid
- self-gravity
  - using Tobias Moldenhauer's correction to Baruteau's implementation
- viscosity
  - physical viscosity
    - explicit/implicit
    - alpha
    - constant
  - numerical visocisty
    - SN
    - TW
- N-body
  - rebound integrator
  - binary
  - initialization with jacobi elements
- frame of reference
  - corotation
  - indirect term calculation
- gravity
  - coupling between N-body and hydro simulation
  - gravity from planet onto gas with potential
  - gravity from gas onto planet with force
  - smoothing
    - in force
    - in potential
    - for all bodies (including star)
    - thickness smoothing
      - cell location
      - planet location for testing
- particles
  - integrator
  - dust
  - dust diffusion
    - random number generation
- opacities
- variable gamma
- radiative processes
  - beta cooling
  - viscous heating
  - irradiative heating
  - radiative cooling
  - radiative transport with flux-limited diffusion
  - QPLUS/QMINUS
- accretion
- aspect ratio
  - isothermal
  - with pressure
- boundary conditions
  - outflow
  - reflective
  - viscous outflow
  - prescribed boundaries
  - mass overflow boundary
- wave damping
- time step calculation


``` {toctree}
:maxdepth: 1
:glob:

*
```