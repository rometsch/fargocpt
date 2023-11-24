# Custom Initial Conditions

This notebook will teach you how to set up a custom gas distribution as initial conditions.

First we create a new directory and change to it.


```python
example_name = "400_custom_initial_conditions"
example_dir = f"example_dirs/{example_name}"
import os
if not os.path.basename(os.getcwd()) == example_name:
    !mkdir -p $example_dir
    os.chdir(example_dir)
repo_root = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
print(f"Current working directory: {os.getcwd()}")
print(f"Repository root directory: {repo_root}")
```

    Current working directory: /home/rometsch/repo/fargocpt/examples/example_dirs/400_custom_initial_conditions
    Repository root directory: /home/rometsch/repo/fargocpt


## Make the code

Make sure the code is built by running make again.

If you have not yet compiled the go, please go to the readme and follow the instructions there.
You can also try to run the following cell directly, but it will only output error messages. This might make debugging harder.


```python
%%timeit -n1 -r1
from sys import platform
if platform in ["linux", "darwin"]:
    !make -j 4 -C $repo_root/src > make.log
else:
    raise RuntimeError(f"Seems like you are not running MacOS or Linux but {platform}. This is unsupported. You are on your own, good luck!")
```

    119 ms ± 0 ns per loop (mean ± std. dev. of 1 run, 1 loop each)


## Preparing a setup file

We'll take the example setup file from the examples directory and modify it in python.
If you want to create setup files for a parameter study, just copy the code and make your own setup creator script.


```python
configfile = "setup.yml"
!cp $repo_root/examples/config.yml $configfile
```

We'll use the `ruamel.yaml` package to read and write the setup file. This can be set up to preserve comments which is very useful if you want to trace your decisions later on.


```python
try:
    import ruamel.yaml
except ImportError:
    raise ImportError("Please install ruamel.yaml with `python3 -m pip install ruamel.yaml`")
yamlloader = ruamel.yaml.YAML()
with open(configfile, "r") as infile:
    config = yamlloader.load(infile)
```


```python
config["MonitorTimestep"] = 0.314 # monitor scalar files around every half orbit
config["Nmonitor"] = 20 # write a snapshot every orbit
config["Nsnapshots"] = 10 # wirte 100 snapshots
# use very low resolution by setting it to 2 cell per scaleheight, cps
config["cps"] = 2

with open(configfile, "w") as outfile:
    yamlloader.dump(config, outfile)
```

## Getting a data skeleton

The next step is to run the code for zero timesteps.
This creates the zeroth and reference snapshots, which we can then modify to our liking.


```python
from fargocpt import run
run(["start", configfile, "-N", "0"], np=2, nt=1, exe=repo_root+"/bin/fargocpt_exe", detach=False)
```

    Running command: mpirun -np 2 --report-pid /tmp/tmp72h2aq3h -x OMP_NUM_THREADS=1 /home/rometsch/repo/fargocpt/bin/fargocpt_exe start setup.yml -N 0
    fargo process pid 1410359
    
    [0] MPI rank #  0 runs as process 1410363
    [0] MPI rank #  0 OpenMP thread #  0 of  1 on cpt-kamino
    [1] MPI rank #  1 runs as process 1410364
    [1] MPI rank #  1 OpenMP thread #  0 of  1 on cpt-kamino
    [0] fargo: This file was compiled on Nov 14 2023, 12:56:40.
    [0] fargo: This version of FARGO used _GNU_SOURCE
    [0] fargo: This version of FARGO used NDEBUG. So no assertion checks!
    [0] Using parameter file setup.yml
    [0] Computing disk quantities within 5.00000e+00 L0 from coordinate center
    [0] BC: Inner composite = reflecting
    [0] BC: Outer composite = reflecting
    [0] BC: Sigma inner = zerogradient
    [0] BC: Sigma outer = zerogradient
    [0] BC: Energy inner = zerogradient
    [0] BC: Energy outer = zerogradient
    [0] BC: Vrad inner = reflecting
    [0] BC: Vrad outer = reflecting
    [0] BC: Vaz inner = keplerian
    [0] BC: Vaz outer = keplerian
    [0] DampingTimeFactor: 1.00000e-01 Outer damping time is computed at radius of 2.50000e+00
    [0] Damping VRadial to reference value at inner boundary.
    [0] Damping VRadial to reference value at outer boundary.
    [0] Damping VAzimuthal to reference value at inner boundary.
    [0] Damping VAzimuthal to reference value at outer boundary.
    [0] Damping SurfaceDensity to reference value at inner boundary.
    [0] Damping SurfaceDensity to reference value at outer boundary.
    [0] Damping Energy to reference value at inner boundary.
    [0] Damping Energy to reference value at outer boundary.
    [0] Radiative diffusion is disabled. Using fixed omega = 1.500000 with a maximum 50000 interations.
    [0] Indirect Term computed as effective Hydro center acceleratrion with shifting the Nbody system to the center.
    [0] Body force on gas computed via potential.
    [0] Using FARGO algorithm for azimuthal advection.
    [0] Using standard forward euler scheme for source terms.
    [0] Cps is set, overwriting Nrad and Naz!
    [0] Grid resolution set using cps = 2.000000
    [0] The grid has (Nrad, Naz) = (74, 251) cells with (1.994113, 1.997395) cps.
    [0] Computing scale height with respect to primary object.
    [0] Using isothermal equation of state. AdiabaticIndex = 1.400.
    [0] Viscosity is of alpha type with alpha = 1.000e-03
    [0] Defaulting to VanLeer flux limiter
    [0] Output information:
    [0]    Output directory: output/out/
    [0]     Number of files: 80
    [0]   Total output size: 0.00 GB
    [0]     Space Available: 31.24 GB
    [0] Initializing 1 RNGs per MPI process.
    [0] Warning : no `radii.dat' file found. Using default.
    [0] The first 1 planets are used to calculate the hydro frame center.
    [0] The mass of the planets used as hydro frame center is 1.000000e+00.
    [0] 2 planet(s) initialized.
    [0] Planet overview:
    [0] 
    [0]  #   | name                    | mass [m0]  | x [l0]     | y [l0]     | vx         | vy         |
    [0] -----+-------------------------+------------+------------+------------+------------+------------+
    [0]    0 | Star                    |          1 |          0 |         -0 |          0 |          0 |
    [0]    1 | Jupiter                 |  0.0009546033 |          1 |          0 |         -0 |   1.000477 |
    [0] 
    [0]  #   | e          | a          | T [t0]     | T [a]      | accreting  | Accretion Type |
    [0] -----+------------+------------+------------+------------+------------+----------------+
    [0]    0 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          0 |   No Accretion |
    [0]    1 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          0 |   No Accretion |
    [0] 
    [0]  #   | Temp [K]   | R [l0]     | irradiates | rampuptime |
    [0] -----+------------+------------+------------+------------+
    [0]    0 |       5778 |  0.0046505 |        yes |          0 |
    [0]    1 |          0 |  4.6505e-05 |         no |          0 |
    [0] 
    [0] Using Tscharnuter-Winkler (1979) artificial viscosity with C = 1.410000.
    [0] Artificial viscosity is used for dissipation.
    [0] Surface density factor: 2.50663
    [0] Tau factor: 0.5
    [0] Tau min: 0.01
    [0] Kappa factor: 1
    [0] Minimum temperature: 2.81162e-05 K = 3.00000e+00
    [0] Maximum temperature: 9.37206e+94 K = 1.00000e+100
    [0] Heating from viscous dissipation is enabled. Using a total factor of 1.
    [0] Cooling (beta) is disabled and reference temperature is floor. Using beta = 10.
    [0] Cooling (radiative) is enabled. Using a total factor of 1.
    [0] S-curve cooling is disabled. 
    [0] CFL parameter: 0.5
    [0] Opacity uses tables from Lin & Papaloizou, 1985
    [0] Particles are disabled.
    [0] Initializing Sigma(r) = 2.25093e-05 = 200 g cm^-2 * [r/(1 AU)]^(-0.5)
    [0] Total disk is mass is 0.000348841 = 6.9366e+29 g.
    [0] Writing output output/out/snapshots/0, Snapshot Number 0, Time 0.000000.
    [0] Writing output output/out/snapshots/reference, Snapshot Number 0, Time 0.000000.
    [0] -- Final: Total Hydrosteps 0, Time 0.00, Walltime 0.01 seconds, Time per Step: 0.00 milliseconds





    0




```python
# %matplotlib widget
from fargocpt import Overview
Overview("output/out",
        vars=["2:Sigma",
                "1:Sigma",
                "1:vazi:rel",
                "1:vrad"]).create()
```


    
![png](400_Custom_Initial_Conditions_files/400_Custom_Initial_Conditions_13_0.png)
    


## Defining helper functions and classes


```python
import numpy as np
from dataclasses import dataclass
from types import SimpleNamespace
import astropy.units as u
import astropy.constants as const
@dataclass
class FargoCPTField:
    outputdir: str
    snapshotid: str
    name: str

    def __post_init__(self):
        self.filename = f"{self.outputdir}/snapshots/{self.snapshotid}/{self.name}.dat"
        self.grid = get_fargo_grid(self.outputdir)
        if self.name == "vrad":
            self.R, self.Phi = np.meshgrid(self.grid.ri, self.grid.phic, indexing="ij")
        elif self.name == "vtheta":
            self.R, self.Phi = np.meshgrid(self.grid.rc, self.grid.phii[:-1], indexing="ij")
        else:
            self.R, self.Phi = np.meshgrid(self.grid.rc, self.grid.phic, indexing="ij")
        self._load()

    def _load(self):
        self._data = np.fromfile(self.filename, dtype=np.float64)
        self._data = self._data.reshape(self.R.shape[0], self.R.shape[1])

    def save(self, altid=None):
        if altid is not None:
            filename = f"{self.outputdir}/snapshots/{altid}/{self.name}.dat"
        self._data.tofile(self.filename)

    @property
    def array(self):
        return self._data

    @array.setter
    def array(self, data):
        if not self._data.shape == data.shape:
            raise ValueError("Shape of data does not match shape of field")
        self._data = data


def get_fargo_grid(outputdir):

    Nrad, Naz = np.genfromtxt(f"{outputdir}/dimensions.dat", usecols=(4,5), dtype=int, unpack=True)

    ri = np.genfromtxt(f"{outputdir}/used_rad.dat")
    phii = np.linspace(0, 2*np.pi, Naz+1)
    Ri, Phii = np.meshgrid(ri, phii, indexing="ij")
    Xi = Ri*np.cos(Phii)
    Yi = Ri*np.sin(Phii)

    rc = 2/3*(ri[1:]**2/(ri[1:]+ri[:-1]) + ri[:-1]) # approx center in polar coords
    phic = 0.5*(phii[1:]+phii[:-1])
    Rc, Phic = np.meshgrid(rc, phic, indexing="ij")
    Xc = Rc*np.cos(Phic)
    Yc = Rc*np.sin(Phic)

    dphi = phii[1] - phii[0]
    dr = ri[1:] - ri[:-1]
    A = 0.5*(Ri[1:,1:]**2 - Ri[:-1,1:]**2)*dphi

    return SimpleNamespace(
        Nrad=Nrad, Naz=Naz,
        ri=ri, phii=phii, Ri=Ri, Phii=Phii, Xi=Xi, Yi=Yi,
        rc=rc, phic=phic, Rc=Rc, Phic=Phic, Xc=Xc, Yc=Yc, 
        dphi=dphi, dr=dr, A=A
    )

def get_fargo_powerlaw_disk(config: dict):

    l0 = u.Unit(config["l0"])
    m0 = u.Unit(config["m0"])
    t0 = u.Unit(config["t0"])
    Temp0 = u.Unit(config["Temp0"])

    Mstar = u.Quantity(config["nbody"][0]["mass"]).to(m0)
    GM = (const.G * Mstar).to(l0**3/t0**2).value
    Rgas = (const.k_B / const.m_p).to(l0**2/t0**2/Temp0).value
    mu = config["mu"]
    r0 = 1
    Sigma0 = u.Quantity(config["Sigma0"])
    if not Sigma0.unit.is_unity():
        Sigma0 = Sigma0.to_value(m0/l0**2)
    betaSigma = -float(config["SigmaSlope"])
    if config["EquationOfState"][:3] == "iso":
        gamma = 1
    else:
        gamma = float(config["AdiabaticIndex"])
    h0 = float(config["AspectRatio"])
    betah = float(config["FlaringIndex"])

    return PowerlawDisk(GM, Rgas, mu, r0, Sigma0, betaSigma, gamma, h0=h0, betah=betah)

def get_fargo_units(outdir: str) -> dict:
    with open(f"{outdir}/units.yml", "r") as f:
        units_info = yaml.safe_load(f)
    units = {}
    for key, val in units_info.items():
        units[key] = u.Unit(f"{val['cgs value']} {val['cgs symbol']}")
    return units

def get_fargo_powerlaw_disk_output(outdir: str):

    config = get_fargo_config(f"{outdir}/snapshots/0/config.yml")
    units = get_fargo_units(outdir)

    l0 = units["length"]
    m0 = units["mass"]
    t0 = units["time"]
    Temp0 = units["temperature"]

    Mstar = u.Quantity(config["nbody"][0]["mass"]).to(m0)
    GM = (const.G * Mstar).to(l0**3/t0**2).value
    Rgas = (const.k_B / const.m_p).to(l0**2/t0**2/Temp0).value
    mu = config["mu"]
    r0 = 1
    Sigma0 = u.Quantity(config["Sigma0"])
    if not Sigma0.unit.is_unity():
        Sigma0 = Sigma0.to_value(m0/l0**2)
    betaSigma = -float(config["SigmaSlope"])
    if config["EquationOfState"][:3] == "iso":
        gamma = 1
    else:
        gamma = float(config["AdiabaticIndex"])
    h0 = float(config["AspectRatio"])
    betah = float(config["FlaringIndex"])

    return PowerlawDisk(GM, Rgas, mu, r0, Sigma0, betaSigma, gamma, h0=h0, betah=betah)

import yaml
def get_fargo_config(filename):
    with open(filename, "r") as f:
        config = yaml.safe_load(f)
    return config


@dataclass
class PowerlawDisk:
    """
    Powerlaw disk model

    Scaleheight is defined with the isothermal sound speed, so H = c_s/sqrt(gamma) / Omega_K

    Parameters
    ----------
    GM : float
        Product of the gravitational constant and the central mass
    Rgas : float
        Gas constant to link c_s and T through the ideal gas law. c_s = sqrt(gamma Rg/mu T). It's Rgas = kB/mH. kB = Boltzman constant. mH = mass of a hydrogen atom.
    mu : float
        Mean molecular weight
    r0 : float
        Reference radius
    Sigma0 : float
        Surface density at r0
    betaSigma : float
        Surface density powerlaw exponent
    T0 : float
        Temperature at r0
    betaT : float
        Temperature powerlaw exponent
    gamma : float
        Adiabatic index
    betaOmega : float (optional)
        Disk angular velocity powerlaw . Defaults to -3/2 for a Keplerian disk.
    """
    GM: float
    Rgas: float
    mu: float
    r0: float
    Sigma0: float
    betaSigma: float
    gamma: float
    T0: float = None
    betaT: float = None
    h0: float = None
    betah: float = None
    betaOmega: float = -3/2

    def __post_init__(self):
        """
        Set derived parameters
        """
        self.OmegaK0 = np.sqrt(self.GM/self.r0**3)
        is_T_set = self.T0 is not None
        is_h_set = self.h0 is not None
        if not (is_T_set or is_h_set):
            raise ValueError("Either T0 or h0 must be set.")
        elif is_T_set and is_h_set:
            raise ValueError("Either T0 or h0 must be set, not both.")
        elif is_T_set:
            self.cs0 = np.sqrt(self.gamma*self.Rgas/self.mu*self.T0)
            self.h0 = self.cs0 / np.sqrt(self.gamma) / self.OmegaK0
        else:
            self.cs0 = self.h0 * np.sqrt(self.gamma) * self.OmegaK0
            self.T0 = self.mu/self.Rgas * self.cs0**2 / self.gamma

        is_betaT_set = self.betaT is not None
        is_betah_set = self.betah is not None
        if not (is_betaT_set or is_betah_set):
            raise ValueError("Either betaT or betah must be set.")
        elif is_betaT_set and is_betah_set:
            raise ValueError("Either betaT or betah must be set, not both.")
        elif is_betaT_set:
            self.betacs = self.betaT/2
            self.betah = self.betaT/2 + 1/2
        else:
            self.betaT = 2*self.betah - 1
            self.betacs = self.betaT/2

        self.betaP = self.betaSigma + self.betaT
        self.P0 = self.Sigma0 * self.cs0**2

    def cs(self, r: np.ndarray):
        return self.cs0 * (r/self.r0)**self.betacs

    def h(self, r: np.ndarray):
        return self.h0 * (r/self.r0)**self.betah

    def Sigma(self, r):
        return self.Sigma0 * (r/self.r0)**self.betaSigma

    def T(self, r: np.ndarray):
        return self.T0 * (r/self.r0)**self.betaT
    
    def P(self, r: np.ndarray):
        return self.P0 * (r/self.r0)**self.betaP

    def OmegaK(self, r: np.ndarray):
        return self.OmegaK0 * (r/self.r0)**self.betaOmega

    def Omega(self, r: np.ndarray):
        """
        Compute the background disk angular velocity profile.

        Use a modification of SB15 Eq. (12) which uses the disk aspect ratio h.
        Omega_0(r) = Omega_K(r) * sqrt{1 + beta_p * h^2}.

        Parameters
        ----------
        r : np.ndarray
            Radial positions

        Returns
        -------
        np.ndarray
        """
        return self.OmegaK(r) * np.sqrt(1 + self.betaP * self.h(r)**2)


    def plot(self, r: np.ndarray):
        """
        Plot an overview of the disks powerlaws.
        """
        mosaic_defs = [
            ["Sigma", "P", "T"],
            ["cs", "h", "Omega"]
        ]
        import matplotlib.pyplot as plt
        fig, axd = plt.subplot_mosaic(mosaic_defs, dpi=150, figsize=(12, 8))

        axd["Sigma"].plot(r, self.Sigma(r))
        axd["Sigma"].set_title("Sigma")
        axd["Sigma"].set_xlabel("r")
        axd["Sigma"].set_ylabel("Sigma")

        axd["P"].plot(r, self.P(r))
        axd["P"].set_title("P")
        axd["P"].set_xlabel("r")
        axd["P"].set_ylabel("P")

        axd["T"].plot(r, self.T(r))
        axd["T"].set_title("T")
        axd["T"].set_xlabel("r")
        axd["T"].set_ylabel("T")
        axd["T"].set_xscale("log")
        axd["T"].set_yscale("log")

        axd["cs"].plot(r, self.cs(r))
        axd["cs"].set_title("cs")
        axd["cs"].set_xlabel("r")
        axd["cs"].set_ylabel("cs")


        axd["h"].plot(r, self.h(r))
        axd["h"].set_title("h")
        axd["h"].set_xlabel("r")
        axd["h"].set_ylabel("h")

        axd["Omega"].plot(r, self.Omega(r))
        axd["Omega"].set_title("Omega")
        axd["Omega"].set_xlabel("r")
        axd["Omega"].set_ylabel("Omega")

        for key, ax in axd.items():
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(which="minor", alpha=0.5,)

        fig.tight_layout()

        return fig, axd

```

## Modify the initial conditions

To modify the data, we first define a helper class to handle the data files.


```python
config = get_fargo_config(configfile)
outdir = "output/out"
disk = get_fargo_powerlaw_disk_output(outdir)
```

### Adding a gap

Let's add a gap by reducing the density and adjusting the azimuthal velocity.


```python
gap_depth = 0.9
gap_width = 0.1
gap_center = 1.0
```


```python
field = FargoCPTField("output/out", "0", "Sigma")
R = field.R
old_density = field.array

f = -gap_depth * np.exp(-0.5*(R-gap_center)**2/gap_width**2)
field.array = old_density * (1+f)

field.save()
```


```python
field = FargoCPTField("output/out", "0", "vazi")
R = field.R
# need to recompute because R is different for vrad
f = -gap_depth * np.exp(-0.5*(R-gap_center)**2/gap_width**2)
# field.array = field.array * ( 1 - 2*R*(R-gap_center)/gap_width**2 * f/(1+f) * disk.cs(R)**2/(disk.gamma*(R*disk.Omega(R))**2) )**0.5
field.array = field.array * ( 1 - 2*R*(R-gap_center)/gap_width**2 * f/(1+f) * disk.h(R)**2 )**0.5

field.save()
```


```python
Overview("output/out",
        vars=["2:Sigma",
                "1:Sigma",
                "1:vazi:rel",
                "1:vrad"]).create()
```


    
![png](400_Custom_Initial_Conditions_files/400_Custom_Initial_Conditions_22_0.png)
    


Now that the initial conditions have been adjusted, we can restart the simulation from snapshot number zero, which will pick up the modified conditions.


```python
run(["auto", configfile], np=2, nt=1, exe=repo_root+"/bin/fargocpt_exe", detach=False)
```

    Running command: mpirun -np 2 --report-pid /tmp/tmp5b2djodn -x OMP_NUM_THREADS=1 /home/rometsch/repo/fargocpt/bin/fargocpt_exe auto setup.yml
    fargo process pid 1410404
    
    [0] MPI rank #  0 runs as process 1410408
    [0] MPI rank #  0 OpenMP thread #  0 of  1 on cpt-kamino
    [1] MPI rank #  1 runs as process 1410409
    [1] MPI rank #  1 OpenMP thread #  0 of  1 on cpt-kamino
    [0] fargo: This file was compiled on Nov 14 2023, 12:56:40.
    [0] fargo: This version of FARGO used _GNU_SOURCE
    [0] fargo: This version of FARGO used NDEBUG. So no assertion checks!
    [0] Using parameter file setup.yml
    [0] Computing disk quantities within 5.00000e+00 L0 from coordinate center
    [0] BC: Inner composite = reflecting
    [0] BC: Outer composite = reflecting
    [0] BC: Sigma inner = zerogradient
    [0] BC: Sigma outer = zerogradient
    [0] BC: Energy inner = zerogradient
    [0] BC: Energy outer = zerogradient
    [0] BC: Vrad inner = reflecting
    [0] BC: Vrad outer = reflecting
    [0] BC: Vaz inner = keplerian
    [0] BC: Vaz outer = keplerian
    [0] DampingTimeFactor: 1.00000e-01 Outer damping time is computed at radius of 2.50000e+00
    [0] Damping VRadial to reference value at inner boundary.
    [0] Damping VRadial to reference value at outer boundary.
    [0] Damping VAzimuthal to reference value at inner boundary.
    [0] Damping VAzimuthal to reference value at outer boundary.
    [0] Damping SurfaceDensity to reference value at inner boundary.
    [0] Damping SurfaceDensity to reference value at outer boundary.
    [0] Damping Energy to reference value at inner boundary.
    [0] Damping Energy to reference value at outer boundary.
    [0] Radiative diffusion is disabled. Using fixed omega = 1.500000 with a maximum 50000 interations.
    [0] Indirect Term computed as effective Hydro center acceleratrion with shifting the Nbody system to the center.
    [0] Body force on gas computed via potential.
    [0] Getting output number of snapshot 0
    [0] Using FARGO algorithm for azimuthal advection.
    [0] Using standard forward euler scheme for source terms.
    [0] Cps is set, overwriting Nrad and Naz!
    [0] Grid resolution set using cps = 2.000000
    [0] The grid has (Nrad, Naz) = (74, 251) cells with (1.994113, 1.997395) cps.
    [0] Computing scale height with respect to primary object.
    [0] Using isothermal equation of state. AdiabaticIndex = 1.400.
    [0] Viscosity is of alpha type with alpha = 1.000e-03
    [0] Defaulting to VanLeer flux limiter
    [0] Output information:
    [0]    Output directory: output/out/
    [0]     Number of files: 80
    [0]   Total output size: 0.00 GB
    [0]     Space Available: 31.23 GB
    [0] Initializing 1 RNGs per MPI process.
    [0] Warning : no `radii.dat' file found. Using default.
    [0] The first 1 planets are used to calculate the hydro frame center.
    [0] The mass of the planets used as hydro frame center is 1.000000e+00.
    [0] 2 planet(s) initialized.
    [0] Planet overview:
    [0] 
    [0]  #   | name                    | mass [m0]  | x [l0]     | y [l0]     | vx         | vy         |
    [0] -----+-------------------------+------------+------------+------------+------------+------------+
    [0]    0 | Star                    |          1 |          0 |         -0 |          0 |          0 |
    [0]    1 | Jupiter                 |  0.0009546033 |          1 |          0 |         -0 |   1.000477 |
    [0] 
    [0]  #   | e          | a          | T [t0]     | T [a]      | accreting  | Accretion Type |
    [0] -----+------------+------------+------------+------------+------------+----------------+
    [0]    0 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          0 |   No Accretion |
    [0]    1 |  6.368246e-17 |          1 |   6.280188 |   0.999548 |          0 |   No Accretion |
    [0] 
    [0]  #   | Temp [K]   | R [l0]     | irradiates | rampuptime |
    [0] -----+------------+------------+------------+------------+
    [0]    0 |       5778 |  0.0046505 |        yes |          0 |
    [0]    1 |          0 |  4.6505e-05 |         no |          0 |
    [0] 
    [0] Using Tscharnuter-Winkler (1979) artificial viscosity with C = 1.410000.
    [0] Artificial viscosity is used for dissipation.
    [0] Surface density factor: 2.50663
    [0] Tau factor: 0.5
    [0] Tau min: 0.01
    [0] Kappa factor: 1
    [0] Minimum temperature: 2.81162e-05 K = 3.00000e+00
    [0] Maximum temperature: 9.37206e+94 K = 1.00000e+100
    [0] Heating from viscous dissipation is enabled. Using a total factor of 1.
    [0] Cooling (beta) is disabled and reference temperature is floor. Using beta = 10.
    [0] Cooling (radiative) is enabled. Using a total factor of 1.
    [0] S-curve cooling is disabled. 
    [0] CFL parameter: 0.5
    [0] Opacity uses tables from Lin & Papaloizou, 1985
    [0] Particles are disabled.
    [0] Initializing Sigma(r) = 2.25093e-05 = 200 g cm^-2 * [r/(1 AU)]^(-0.5)
    [0] Total disk is mass is 0.000348841 = 6.9366e+29 g.
    [0] Loading polargrinds for damping...
    [0] Reading file 'output/out/snapshots/reference/Sigma.dat' with 148592 bytes.
    [0] Reading file 'output/out/snapshots/reference/vrad.dat' with 150600 bytes.
    [0] Reading file 'output/out/snapshots/reference/vazi.dat' with 148592 bytes.
    [0] Loading polargrinds at t = 0...
    [0] Reading file 'output/out/snapshots/0/Sigma.dat' with 148592 bytes.
    [0] Reading file 'output/out/snapshots/0/vrad.dat' with 150600 bytes.
    [0] Reading file 'output/out/snapshots/0/vazi.dat' with 148592 bytes.
    [0] Restarting planetary system...
    [0] Loading planets ...[0]  done
    [0] Loading rebound ...[0]  done
    [0] Finished restarting planetary system.
    [0] Writing output output/out/snapshots/1, Snapshot Number 1, Time 6.280000.
    [0] Writing output output/out/snapshots/2, Snapshot Number 2, Time 12.560000.
    [0] Writing output output/out/snapshots/3, Snapshot Number 3, Time 18.840000.
    [0] Writing output output/out/snapshots/4, Snapshot Number 4, Time 25.120000.
    [0] Writing output output/out/snapshots/5, Snapshot Number 5, Time 31.400000.
    [0] Writing output output/out/snapshots/6, Snapshot Number 6, Time 37.680000.
    [0] Writing output output/out/snapshots/7, Snapshot Number 7, Time 43.960000.
    [0] Writing output output/out/snapshots/8, Snapshot Number 8, Time 50.240000.
    [0] Writing output output/out/snapshots/9, Snapshot Number 9, Time 56.520000.
    [0] Writing output output/out/snapshots/10, Snapshot Number 10, Time 62.800000.
    [0] -- Final: Total Hydrosteps 1728, Time 62.80, Walltime 4.01 seconds, Time per Step: 2.32 milliseconds





    0



And finally, let's see how this system evolved.

Uncomment the first line to get an interactive widget, if your jupyter installation supports it.


```python
%matplotlib inline
Overview("output/out",
        vars=["2:Sigma",
                "1:Sigma",
                "1:vazi:rel",
                "1:vrad"]).create()
```


    
![png](400_Custom_Initial_Conditions_files/400_Custom_Initial_Conditions_26_0.png)
    



```python

```
