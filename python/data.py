import yaml
import os

import numpy as np

from os.path import join as joinpath

from dataclasses import dataclass, field

from astropy.units import Unit, Quantity

from typing import List

from functools import lru_cache


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
        """ Return a meshgrid of the radial and azimuthal coordinates where the data is defined.

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


@lru_cache(20)
def _load_text_data_variables(filepath, timestamp):
    # load all variable definitions from a text file
    # which contains the variable names and colums in its header.
    # each variable is indicated by
    # "#variable: {column number} | {variable name} | {unit}"
    found_variables = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line[0] != "#":
                break
            identifier = "#variable:"
            if line[:len(identifier)] == identifier:
                col, name, unitstr = [
                    s.strip() for s in line[len(identifier):].split("|")
                ]
                found_variables[name] = (col, unitstr)
    return found_variables

def load_text_data_variables(filepath):
    return _load_text_data_variables(filepath, os.path.getmtime(filepath))


@lru_cache(20)
def _load_data(filepath, timestamp):
    return np.genfromtxt(filepath).T

def load_data(filepath):
    """ Load data from a text file and chache the results. 
    The chache is deleted once the file is modified."""
    return _load_data(filepath, os.path.getmtime(filepath))


def load_text_data_file(filepath, varname, Nmax=np.inf):
    # get data
    variables = load_text_data_variables(filepath)
    col = int(variables[varname][0])
    unit_str = variables[varname][1]
    unit_str = unit_str.replace("1/s", "s-1")
    unit = Unit(unit_str)
    file_data = load_data(filepath)
    data = file_data[col] * unit
    if data.isscalar:
        data = Quantity([data])
    return data

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
    grid: Grid
    target_units: Units = None

    def __post_init__(self):
        self._load_info()
    
    def _load_info(self):
        info_file = joinpath(self.output_dir, 'info1D.yml')
        try:
            with open(info_file, "r") as f:
                self._info_dict = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError("Could not find 'info1D.yml' in simulatiuon output ('{}')".format(self.output_dir))
        if self._info_dict is None:
            self._info_dict = {}

    @property
    def var_names(self):
        return [k for k in self._info_dict.keys()]
    
    def __repr__(self) -> str:
        rv = "FargoCPT data Vars1D\n"
        rv += "====================\n"
        rv += f"  output_dir: {self.output_dir}\n"
        rv += f"  target_units" + ("= None" if self.target_units is None else ": Units") + "\n"
        rv += f"  grid: Grid\n"
        rv += f"  var_names:\n"
        for var_name in self.var_names:
            rv += f"    {var_name}\n"
        rv += "====================\n"
        return rv
    
    def _return_radius(self, varname):
        info = self._info_dict[varname]
        rad_interface = info["on_radial_interface"]
        if rad_interface:
            r = self.grid.radi
        else:
            r = self.grid.radc
        return r

    def _return_data(self, varname, Nsnapshot, data_slice):

        info = self._info_dict[varname]
        unit = Unit(info['unit'])

        filename = info["filename"]
        datafile = joinpath(self.output_dir, 'snapshots', f'{Nsnapshot}', filename)
        data = np.fromfile(datafile, dtype=np.float64)[data_slice]
        data = data * unit

        return data

    def get(self, varname, Nsnapshot, include_grid=True):
        return self.average(varname, Nsnapshot, include_grid=include_grid)

    def average(self, varname, Nsnapshot, include_grid=True):
        avg_slice = slice(1,None,4)
        data = self._return_data(varname, Nsnapshot, avg_slice)
        if include_grid:
            r = self._return_radius(varname)
            return r, data
        else:
            return data
    
    def min(self, varname, Nsnapshot, include_grid=True):
        max_slice = slice(2,None,4)
        data = self._return_data(varname, Nsnapshot, max_slice)
        if include_grid:
            r = self._return_radius(varname)
            return r, data
        else:
            return data
    
    def max(self, varname, Nsnapshot, include_grid=True):
        min_slice = slice(3,None,4)
        data = self._return_data(varname, Nsnapshot, min_slice)
        if include_grid:
            r = self._return_radius(varname)
            return r, data
        else:
            return data
    
@dataclass
class Vars2D:
    output_dir: str
    grid: Grid
    target_units: Units = None

    def __post_init__(self):
        info_file = joinpath(self.output_dir, 'info2D.yml')
        try:
            with open(info_file, "r") as f:
                self._info_dict = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError("Could not find 'info2D.yml' in simulatiuon output ('{}')".format(self.output_dir))
        if self._info_dict is None:
            self._info_dict = {}

    @property
    def var_names(self):
        return [k for k in self._info_dict.keys()]
    
    def __repr__(self) -> str:
        rv = "FargoCPT data Vars2D\n"
        rv += "====================\n"
        rv += f"  output_dir: {self.output_dir}\n"
        rv += f"  target_units" + ("= None" if self.target_units is None else ": Units") + "\n"
        rv += f"  grid: Grid\n"
        rv += f"  var_names:\n"
        for var_name in self.var_names:
            rv += f"    {var_name}\n"
        rv += "====================\n"
        return rv

    def meshgrid(self, varname, Nsnapshot):
        info = self._info_dict[varname]
        rad_interface = info["on_radial_interface"]
        azi_interface = info["on_azimuthal_interface"]
        return self.grid.meshgrid(intr=rad_interface, intf=azi_interface)

    def meshgrid_plot(self, varname, Nsnapshot):
        info = self._info_dict[varname]
        rad_interface = info["on_radial_interface"]
        azi_interface = info["on_azimuthal_interface"]
        return self.grid.meshgrid_plot(intr=rad_interface, intf=azi_interface)

    def get(self, varname, Nsnapshot, include_grid=True, grid_for_plot=False):
        if not varname in self.var_names:
            raise KeyError(f"Unknown variable '{varname}'")

        info = self._info_dict[varname]
        unit = Unit(info['unit'])

        Nr = info["Nrad"]
        Naz = info["Nazi"]
        filepath = joinpath(self.output_dir, "snapshots", f"{Nsnapshot}", info["filename"])
        rv = np.fromfile(filepath).reshape(Nr, Naz) * unit

        if include_grid:
            if grid_for_plot:
                r, phi = self.meshgrid_plot(varname, Nsnapshot)
            else:
                r, phi = self.meshgrid(varname, Nsnapshot)
            return r, phi, rv
        else:
            return rv
        
    def _return_with_radius(self, varname, data, include_grid):
        if include_grid:
            if self._info_dict[varname]["on_radial_interface"]:
                r = self.grid.radi
            else:
                r = self.grid.radc
            return r, data
        else:
            return data

    def average(self, varname, Nsnapshot, include_grid=True):
        data = np.average(self.get(varname, Nsnapshot, include_grid=False), axis=1)
        return self._return_with_radius(varname, data, include_grid)

    def min(self, varname, Nsnapshot, include_grid=True):
        data = np.min(self.get(varname, Nsnapshot, include_grid=False), axis=1)
        return self._return_with_radius(varname, data, include_grid)
        
    def max(self, varname, Nsnapshot, include_grid=True):
        data = np.max(self.get(varname, Nsnapshot, include_grid=False), axis=1)
        return self._return_with_radius(varname, data, include_grid)
        


@dataclass
class Particles:
    output_dir: str
    target_units: Units = None
    _last_snapshot = None
    _last_data = None

    def __post_init__(self):
        self._load_info()

    def _load_info(self):
        info_file = joinpath(self.output_dir, 'infoParticles.yml')
        try:
            with open(info_file, "r") as f:
                self._info_dict = yaml.safe_load(f)
        except FileNotFoundError:
            raise FileNotFoundError("Could not find 'infoParticles.yml' in simulatiuon output ('{}')".format(self.output_dir))
    
        self.var_names = [k for k in self._info_dict['variables'].keys()]
        self._dtypes = []
        for key, info in self._info_dict['variables'].items():
            dtype = info['type']
            if dtype == 'double':
                self._dtypes.append((key, np.dtype(np.dtype('f8'))))
            elif dtype == 'unsigned long':
                self._dtypes.append((key, 'u8'))
            else:
                raise ValueError(f'Unknown type {dtype}')
            
        if "r" in self.var_names and "phi" in self.var_names:
            self.var_names.append("x")
            self.var_names.append("y")
        elif "x" in self.var_names and "y" in self.var_names:
            self.var_names.append("r")
            self.var_names.append("phi")

    def _load_snapshot(self, Nsnapshot):
        # cache the last read snapshot
        filename = joinpath(self.output_dir, "snapshots", f"{Nsnapshot}", "particles.dat")
        if self._last_snapshot == Nsnapshot:
            res = self._last_data
        else:
            res = np.fromfile(filename, dtype=self._dtypes)
            self._last_snapshot = Nsnapshot
            self._last_data = res
        return res

    def _handle_xy_rphi(self, varname, res):
        r_direct = varname == 'r' and 'r' in self._info_dict['variables']
        phi_direct = varname == 'phi' and 'phi' in self._info_dict['variables']
        x_direct = varname == 'x' and 'x' in self._info_dict['variables']
        y_direct = varname == 'y' and 'y' in self._info_dict['variables']
        if r_direct or phi_direct or x_direct or y_direct:
            unit = Unit(self._info_dict['variables'][varname]['unit'])
            rv = res[varname]*unit
            return rv

        # otherwise calculate the values
        if varname == 'r':
            unit = Unit(self._info_dict['variables']["x"]['unit'])
            rv = np.sqrt(res['x']**2 + res['y']**2)*unit
        elif varname == 'phi':
            rv = np.arctan2(res['y'], res['x'])
        elif varname == 'x':
            unit = Unit(self._info_dict['variables']["r"]['unit'])
            rv = res['r']*np.cos(res['phi'])*unit
        elif varname == 'y':
            unit = Unit(self._info_dict['variables']["r"]['unit'])
            rv = res['r']*np.sin(res['phi'])*unit
        return rv
    

    def get(self, varname, Nsnapshot):
        if not varname in self.var_names:
            raise KeyError(f"Unknown particle variable '{varname}'")

        res = self._load_snapshot(Nsnapshot)

        if varname in ['x', 'y', 'r', 'phi']:
            rv = self._handle_xy_rphi(varname, res)
        else:
            unit = Unit(self._info_dict['variables'][varname]['unit'])
            rv = res[varname]*unit

        return rv

    def __repr__(self) -> str:
        rv = "FargoCPT data Particles\n"
        rv += "====================\n"
        rv += f"  output_dir: {self.output_dir}\n"
        rv += f"  target_units" + ("= None" if self.target_units is None else ": Units") + "\n"
        rv += f"  var_names:\n"
        for var_name in self.var_names:
            rv += f"    {var_name}\n"
        rv += "====================\n"
        return rv

@dataclass
class Hydro:
    output_dir: str
    units: Units
    target_units: Units = None
    grid: Grid = None
    timestepping: Timestepping = None
    scalars: Scalars = None
    vars1D: Vars1D = None

    def __post_init__(self):
        self.timestepping = Timestepping(0, joinpath(self.output_dir, 'monitor', 'timestepLogging.dat'))
        self.scalars = Quantities(0, joinpath(self.output_dir, 'monitor', 'Quantities.dat'))
        self.grid = Grid()
        self.grid.load_files(self.output_dir, self.units, target_units=self.target_units)
        self.vars1D = Vars1D(self.output_dir, self.grid, target_units=self.target_units)
        self.vars2D = Vars2D(self.output_dir, self.grid, target_units=self.target_units)

    def __repr__(self) -> str:
        rv = "FargoCPT data Hydro\n"
        rv += "====================\n"
        rv += f"  output_dir: {self.output_dir}\n"
        rv += f"  units: Units\n"
        rv += f"  target_units" + ("= None" if self.target_units is None else ": Units") + "\n"
        rv += f"  grid: Grid\n"
        rv += f"  timestepping: Scalar\n"
        rv += f"  scalars: Scalar\n"
        rv += f"  vars1D: Vars1D\n"
        rv += f"  vars2D: Vars2D\n"
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
    particles: Particles = None

    def __post_init__(self):
        if not os.path.exists(self.output_dir):
            raise FileNotFoundError(f"Could not find output directory '{self.output_dir}'")

        if self.units is None:
            self.units = Units()
            self.units.load_file(joinpath(self.output_dir , 'units.yml'))

        self._load_snapshots()

        self.hydro = Hydro(self.output_dir, self.units, target_units=self.target_units)
        self._load_nbody()

        self.params = Params(self.output_dir)
        self._load_particles()

    def _load_nbody(self):
        for n in range(1,100):
            path = joinpath(self.output_dir, 'monitor', f'planet{n}.dat')
            if os.path.exists(path):
                nbody = Nbody(n, filepath=path)
                self.nbody.append(nbody)
            else:
                break

    def _load_particles(self):
        info_file = joinpath(self.output_dir, 'infoParticles.yml')
        if os.path.exists(info_file):
            self.particles = Particles(self.output_dir, target_units=self.target_units)

    def _load_snapshots(self):
        filename = joinpath(self.output_dir, 'snapshots', 'list.txt')
        self.snapshots = []
        self.special_snapshots = []
        with open(filename, 'r') as f:
            for line in f:
                try:
                    self.snapshots.append(int(line))
                except ValueError:
                    self.special_snapshots.append(line.strip())
        if os.path.exists(joinpath(self.output_dir, 'snapshots', 'reference')):
            self.special_snapshots.append('reference')

        self.snapshot_time = load_text_data_file(joinpath(self.output_dir, 'snapshots', 'timeSnapshot.dat'), 'time')
        self.monitor_number = [int(n) for n in load_text_data_file(joinpath(self.output_dir, 'snapshots', 'timeSnapshot.dat'), 'monitor number')]

    def snapshot_time(self, Nsnapshot):
        filename = joinpath(self.output_dir, 'snapshots', f'{Nsnapshot}', 'time.dat')
        return load_text_data_file(filename, 'time')


    def __repr__(self) -> str:
        rv = "FargoCPT data Loader\n"
        rv += "====================\n"
        rv += f"  output_dir: {self.output_dir}\n"
        rv += f"  snapshots: {self.snapshots[0]} ... {self.snapshots[-1]}\n"
        if len(self.special_snapshots) > 0:
            rv += f"  special_snapshots: {self.special_snapshots}\n"
        rv += f"  snapshot_time: {self.snapshot_time[0]} ... {self.snapshot_time[-1]}\n"
        rv += f"  monitor_number: {self.monitor_number[0]} ... {self.monitor_number[-1]}\n"
        rv += f"  units: Units\n"
        rv += f"  target_units" + ("= None" if self.target_units is None else ": Units") + "\n"
        rv += f"  hydro: Hydro\n"
        rv += f"  nbody: Nbody\n"
        rv += f"  params: Params\n"
        rv += f"  particles" + (" = None" if self.particles is None else ": Particles") + "\n"
        rv += "====================\n"
        return rv