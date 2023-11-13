import yaml
import os

import numpy as np

from os.path import join as joinpath

from dataclasses import dataclass, field

from astropy.units import Unit

from load_scalar import load_text_data_file, load_text_data_variables

from typing import List


@dataclass
class Units:
    length: Unit = None
    time: Unit = None
    mass: Unit = None
    temperature: Unit = None

    def load_file(self, path):
        with open(path, 'r') as f:
            unitdata = yaml.safe_load(f)
        for key, value in unitdata.items():
            setattr(self, key, Unit(value["unit"]))

@dataclass
class Grid:
    radi = None
    phii = None
    rmin = None
    rmax = None
    Naz = None
    Nrad = None
    spacing = None
    
    @property
    def radc(self):
        if self.radi is None:
            return None
        return 0.5*(self.radi[1:] + self.radi[:-1])
        
    @property
    def phic(self):
        if self.phii is None:
            return None
        return 0.5*(self.phii[1:] + self.phii[:-1])
    
    @property
    def dphi(self):
        if self.phii is None:
            return None
        return self.phii[1] - self.phii[0]
    
    @property
    def drad(self):
        if self.radi is None:
            return None
        return self.radi[1:] - self.radi[:-1]

    def meshgrid_plot(self, intr = False, intf = False):
        """ Return a meshgrid of the radial and azimuthal coordinates for plotting purposes.

        This function is similar to `meshgrid`, but it extents the radial coordinates by one in case if intr is True.
        In case of intf is True, the cell centers are used with the first value repeated + 2pi.
        The result can be used to plot the data on a pcolormesh plot, which expects one more cell than the data has.

        Parameters
        ----------
        intr : bool
            If True, the returned radial coordinates will be the radial cell interfaces.
        intf : bool
            If True, the returned azimuthal coordinates will be the azimuthal cell interfaces.
        
        Returns
        -------
        r, phi : tuple of np.ndarray
            The radial and azimuthal coordinates.
        """
        if intr:
            r = np.append([self.radc[0]-self.drad[0]], self.radc)
            r = np.append(r, [self.radc[-1]+self.drad[-1]])
        else:
            r = self.radi
        if intf:
            phi = np.append(self.phic, [self.phic[0]+2*np.pi])
        else:
            phi = self.phii
        return np.meshgrid(r, phi, indexing='ij')

    def meshgrid(self, intr = False, intf = False):
        """ Return a meshgrid of the radial and azimuthal coordinates.

        Parameters
        ----------
        intr : bool
            If True, the returned radial coordinates will be the radial cell interfaces.
        intf : bool
            If True, the returned azimuthal coordinates will be the azimuthal cell interfaces.
        
        Returns
        -------
        r, phi : tuple of np.ndarray
            The radial and azimuthal coordinates.
        """
        if intr:
            r = self.radi
        else:
            r = self.radc
        if intf:
            phi = self.phii
        else:
            phi = self.phic
        return np.meshgrid(r, phi, indexing='ij')

    def load_files(self, output_dir, units, target_units=None):
        dimensions_file = joinpath(output_dir, 'dimensions.dat')
        with open(dimensions_file, 'r') as file:
            next(file)  # Skip the header line
            data_line = next(file)  # Read the data line
            fields = data_line.split()  # Split the data line into fields

        nrad, naz, nghrad, nghaz = map(int, fields[4:8])
        phimin = 0
        phimax = 2*np.pi
        self.phii = np.linspace(phimin, phimax, naz+1, endpoint=True)
        # Convert the fields to the appropriate types

        self.Nrad = nrad
        self.Naz = naz

        radius_file = joinpath(output_dir, 'used_rad.dat')
        self.radi = np.loadtxt(radius_file)*units.length
        if target_units is not None:
            self.radi = self.radi.to(target_units.length)

        self.rmin = self.radi[0]
        self.rmax = self.radi[-1]
        self.spacing = fields[8]

@dataclass
class Nbody:
    id: int
    filepath: str
    _target_units: Units = None
    _varnames = ["time", "snapshot_number", "monitor_number", "x", "y", "vx", "vy", "mass", "physical_time", "omega_frame", "mdcp", "eccentricity", "angular_momentum", "semi_major_axis", "omega_kepler", "mean_anomaly", "eccentric_anomaly", "true_anomaly", "pericenter_angle", "torque", "accreted_mass", "accretion_rate"]
    
    def _return_var(self, varname):
        rv = load_text_data_file(self.filepath, varname)
        return rv

    @property
    def time(self):
        return self._return_var('time')
    
    @property
    def snapshot_number(self):
        return self._return_var('snapshot number')

    @property
    def monitor_number(self):
        return self._return_var('monitor number')

    @property
    def x(self):
        return self._return_var('x')

    @property
    def y(self):
        return self._return_var('y')

    @property
    def vx(self):
        return self._return_var('vx')

    @property
    def vy(self):
        return self._return_var('vy')

    @property
    def mass(self):
        return self._return_var('mass')

    @property
    def physical_time(self):
        return self._return_var('physical time')

    @property
    def omega_frame(self):
        return self._return_var('omega frame')

    @property
    def mdcp(self):
        return self._return_var('mdcp')

    @property
    def eccentricity(self):
        return self._return_var('eccentricity')

    @property
    def angular_momentum(self):
        return self._return_var('angular momentum')

    @property
    def semi_major_axis(self):
        return self._return_var('semi-major axis')

    @property
    def omega_kepler(self):
        return self._return_var('omega kepler')

    @property
    def mean_anomaly(self):
        return self._return_var('mean anomaly')

    @property
    def eccentric_anomaly(self):
        return self._return_var('eccentric anomaly')

    @property
    def true_anomaly(self):
        return self._return_var('true anomaly')

    @property
    def pericenter_angle(self):
        return self._return_var('pericenter angle')

    @property
    def torque(self):
        return self._return_var('torque')

    @property
    def accreted_mass(self):
        return self._return_var('accreted mass')

    @property
    def accretion_rate(self):
        return self._return_var('accretion rate')


@dataclass
class Scalars:
    id: int
    filepath: str
    _target_units: Units = None
    _varnames = None

    def _return_var(self, varname):
        rv = load_text_data_file(self.filepath, varname)
        return rv

    def __repr__(self):
        rv = "FargoCPT data Scalars\n"
        rv += "====================\n"
        rv += f"  filepath: {self.filepath}\n"
        rv += f"  varnames:\n"
        for varname in self._varnames:
            rv += f"    {varname.replace(' ', '_')}\n"
        rv += "====================\n"
        return rv

class Timestepping(Scalars):
    _varnames = ["snapshot_number", "monitor_number", "hydrostep_number", "Number_of_Hydrosteps_in_last_monitor_timestep", "time", "walltime", "walltime_per_hydrostep", "mean_dt", "min_dt", "std_dev_dt"]
    
    @property
    def snapshot_number(self):
        return self._return_var('snapshot number')

    @property
    def monitor_number(self):
        return self._return_var('monitor number')

    @property
    def hydrostep_number(self):
        return self._return_var('hydrostep number')

    @property
    def Number_of_Hydrosteps_in_last_monitor_timestep(self):
        return self._return_var('Number of Hydrosteps in last monitor timestep')

    @property
    def time(self):
        return self._return_var('time')

    @property
    def walltime(self):
        return self._return_var('walltime')

    @property
    def walltime_per_hydrostep(self):
        return self._return_var('walltime per hydrostep')

    @property
    def mean_dt(self):
        return self._return_var('mean dt')

    @property
    def min_dt(self):
        return self._return_var('min dt')

    @property
    def std_dev_dt(self):
        return self._return_var('std dev dt')


class Quantities(Scalars):
    _varnames = ["snapshot number", "monitor number", "time", "mass", "radius", "angular momentum", "total energy", "internal energy", "kinematic energy", "potential energy", "radial kinetic energy", "azimuthal kinetic energy", "eccentricity", "periastron", "viscous dissipation", "luminosity", "pdivv", "inner boundary mass inflow", "inner boundary mass outflow", "outer boundary mass inflow", "outer boundary mass outflow", "wave damping inner mass creation", "wave damping inner mass removal", "wave damping outer mass creation", "wave damping outer mass removal", "density floor mass creation", "aspect", "indirect term nbody x", "indirect term nbody y", "indirect term disk x", "indirect term disk y", "frame angle", "advection torque", "viscous torque", "gravitational torque"]
    
    @property
    def snapshot_number(self):
        return self._return_var('snapshot number')

    @property
    def monitor_number(self):
        return self._return_var('monitor number')

    @property
    def time(self):
        return self._return_var('time')

    @property
    def mass(self):
        return self._return_var('mass')

    @property
    def radius(self):
        return self._return_var('radius')

    @property
    def angular_momentum(self):
        return self._return_var('angular momentum')

    @property
    def total_energy(self):
        return self._return_var('total energy')

    @property
    def internal_energy(self):
        return self._return_var('internal energy')

    @property
    def kinematic_energy(self):
        return self._return_var('kinematic energy')

    @property
    def potential_energy(self):
        return self._return_var('potential energy')

    @property
    def radial_kinetic_energy(self):
        return self._return_var('radial kinetic energy')

    @property
    def azimuthal_kinetic_energy(self):
        return self._return_var('azimuthal kinetic energy')

    @property
    def eccentricity(self):
        return self._return_var('eccentricity')

    @property
    def periastron(self):
        return self._return_var('periastron')

    @property
    def viscous_dissipation(self):
        return self._return_var('viscous dissipation')

    @property
    def luminosity(self):
        return self._return_var('luminosity')

    @property
    def pdivv(self):
        return self._return_var('pdivv')

    @property
    def inner_boundary_mass_inflow(self):
        return self._return_var('inner boundary mass inflow')

    @property
    def inner_boundary_mass_outflow(self):
        return self._return_var('inner boundary mass outflow')

    @property
    def outer_boundary_mass_inflow(self):
        return self._return_var('outer boundary mass inflow')

    @property
    def outer_boundary_mass_outflow(self):
        return self._return_var('outer boundary mass outflow')

    @property
    def wave_damping_inner_mass_creation(self):
        return self._return_var('wave damping inner mass creation')

    @property
    def wave_damping_outer_mass_removal(self):
        return self._return_var('wave damping outer mass removal')

    @property
    def density_floor_mass_creation(self):
        return self._return_var('density floor mass creation')

    @property
    def aspect(self):
        return self._return_var('aspect')

    @property
    def indirect_term_nbody_x(self):
        return self._return_var('indirect term nbody x')

    @property
    def indirect_term_nbody_y(self):
        return self._return_var('indirect term nbody y')

    @property
    def indirect_term_disk_x(self):
        return self._return_var('indirect term disk x')

    @property
    def indirect_term_disk_y(self):
        return self._return_var('indirect term disk y')

    @property
    def frame_angle(self):
        return self._return_var('frame angle')

    @property
    def advection_torque(self):
        return self._return_var('advection torque')

    @property
    def viscous_torque(self):
        return self._return_var('viscous torque')

    @property
    def gravitational_torque(self):
        return self._return_var('gravitational torque')

@dataclass
class Vars1D:
    output_dir: str
    units: Units
    target_units: Units = None

    def __post_init__(self):
        self._load_info()
        self._parse_info()
    
    def _load_info(self):
        info_file = joinpath(self.output_dir, 'info1D.yml')
        try:
            with open(self.datadir_path(info_file), "r") as f:
                self._info_dict = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError("Could not find 'info2D.yml' in simulatiuon output ('{}')".format(self.output_dir))

    def _parse_info(self):
        pass

@dataclass
class Hydro:
    output_dir: str
    units: Units
    target_units: Units = None
    grid: Grid = None
    timestepping: Timestepping = None
    scalars: Scalars = None

    def __post_init__(self):
        self.timestepping = Timestepping(0, joinpath(self.output_dir, 'monitor', 'timestepLogging.dat'))
        self.scalars = Quantities(0, joinpath(self.output_dir, 'monitor', 'Quantities.dat'))
        self.grid = Grid()
        self.grid.load_files(self.output_dir, self.units, target_units=self.target_units)

    def __repr__(self) -> str:
        rv = "FargoCPT data Hydro\n"
        rv += "====================\n"
        rv += f"  output_dir: {self.output_dir}\n"
        rv += f"  units: Units\n"
        rv += f"  target_units" + ("= None" if self.target_units is None else ": Units") + "\n"
        rv += f"  grid: Grid\n"
        rv += f"  timestepping: Scalar\n"
        rv += f"  scalars: Scalar\n"
        rv += "====================\n"
        return rv

class Params:
    """ Class to load the parameters from the output directory supporting non-case sensitive keys."""
    def __init__(self, output_dir, snapshot=None):
        self.output_dir = output_dir
        self.snapshot = snapshot
        self._load_params()

    def _get_newest_param_file(self):
        directory = joinpath(self.output_dir, 'parameters')
        full_paths = [os.path.join(directory, filename) for filename in os.listdir(directory)]
        full_paths = [path for path in full_paths if os.path.isfile(path)]
        newest_file = max(full_paths, key=os.path.getmtime)
        return newest_file

    def _load_params(self):
        if self.snapshot is None:
            self._param_filepath = self._get_newest_param_file()
        else:
            self._param_filepath = joinpath(self.output_dir, 'snapshots', f'{self.snapshot}', 'config.yml')
        with open(self._param_filepath, 'r') as f:
            self.params = yaml.safe_load(f)

        self.lkeys = {key.lower(): key for key in self.params.keys()}

    def __getitem__(self, key):
        return self.params[self.lkeys[key.lower()]]
    
    def __repr__(self) -> str:
        rv = "FargoCPT data Params\n"
        rv += "====================\n"
        rv += f"  filename: {self._param_filepath}\n"
        rv += f"  params:\n"
        for key, value in self.params.items():
            rv += f"    {key}: {value}\n"
        rv += "====================\n"
        return rv

@dataclass
class Loader:
    output_dir: str
    units: Units = None
    target_units: Units = None
    hydro: Hydro = None
    nbody: List[Nbody] = field(default_factory=list)
    params: Params = None

    def __post_init__(self):
        if self.units is None:
            self.units = Units()
            self.units.load_file(joinpath(self.output_dir , 'units.yml'))

        self.hydro = Hydro(self.output_dir, self.units, target_units=self.target_units)
        self.load_nbody()

        self.params = Params(self.output_dir)

    def load_nbody(self):
        for n in range(1,100):
            path = joinpath(self.output_dir, 'monitor', f'planet{n}.dat')
            if os.path.exists(path):
                nbody = Nbody(n, filepath=path)
                self.nbody.append(nbody)
            else:
                break

    def __repr__(self) -> str:
        rv = "FargoCPT data Loader\n"
        rv += "====================\n"
        rv += f"  output_dir: {self.output_dir}\n"
        rv += f"  units: Units\n"
        rv += f"  target_units" + ("= None" if self.target_units is None else ": Units") + "\n"
        rv += f"  hydro: Hydro\n"
        rv += f"  nbody: Nbody\n"
        rv += f"  params: Params\n"
        rv += f"  particles: Particles\n"
        rv += "====================\n"
        return rv