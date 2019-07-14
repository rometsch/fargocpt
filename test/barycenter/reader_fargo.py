# A library providing functions to read fargo data
import numpy as np
import astropy.units as u
import os

# A dict containing the scalar variables.
# The entries have the following meaning
# file : file to load the data from
# datacol : column in the file where the data is stored
# timecol : column in the file where the time is stored
# unitpowers : powers of the baseunits corresponding to the unit of the observable
scalarVars = {
    'mass' : {
        'file' : 'Quantities.dat',
        'datacol' : 1,
        'timecol' : 0,
        'unitpowers' : {'mass' : 1} },
     'ekin_radial' : {
        'file' : 'Quantities.dat',
        'datacol' : 10,
        'timecol' : 0,
         'unitpowers' : {'mass' : 1, 'length' : 2, 'time' : -2} },
    'ekin_azimuthal' : {
        'file' : 'Quantities.dat',
        'datacol' : 11,
        'timecol' : 0,
        'unitpowers' : {'mass' : 1, 'length' : 2, 'time' : -2} },
    'lostmass_inner' : {
        'file' : 'bigplanet1.dat',
        'datacol' : 6,
        'timecol' : 7,
        'unitpowers' : {'mass' : 1} }

    }

# A dict containing the scalar variables for each planet.
# There is one file for each planet.
# In the following notation, the planets are indexed starting from 1
# and the offset specifies the filename indexing convention.
# E.g. the first planet corresponds to the file 'orbit0.dat'
# and has an offset of -1.
# For the planet files, the file string contains a '{}' to be used together with format(nplanet+offset)
# The other entries are as above.

scalarVarsPlanet = {
    'a' : {
        'file' : 'bigplanet{}.dat',
        'offset' : 0,
        'datacol' : 13,
        'timecol' : 7,
        'unitpowers' : {'length' : 1} },
    'ecc' : {
        'file' : 'bigplanet{}.dat',
        'offset' : 0,
        'datacol' : 11,
        'timecol' : 7,
        'unitpowers' : {} },
    'x' : {
        'file' : 'bigplanet{}.dat',
        'offset' : 0,
        'datacol' : 1,
        'timecol' : 7,
        'unitpowers' : {'length' : 1} },
    'y' : {
        'file' : 'bigplanet{}.dat',
        'offset' : 0,
        'datacol' : 2,
        'timecol' : 7,
        'unitpowers' : {'length' : 1} },
    'mass' : {
        'file' : 'bigplanet{}.dat',
        'offset' : 0,
        'datacol' : 5,
        'timecol' : 7,
        'unitpowers' : {'mass' : 1} }
    }

# A dict containing the 1D output variables.
units1D = {
	'gasdens': u.g/u.cm**2,
	'gasEccentricity': 1,
	'gasTemperature': u.K,
	'gasenergy': u.erg,
	'gasvrad' : u.cm/u.s,
	'gasvtheta' : u.cm/u.s,
	'gasTemperature' : u.K
    }

# Aliases to match variable names in the output data with previously used names
# or names from other codes
aliases = {
    'ekin_radial' : "radial kinetic energy",
    'ekin_azimuthal' : "azimuthal kinetic energy",
    "a" : "semi-major axis",
    "ecc" : "eccentricity"
}


def loadCoarseOutputTimes(dataDir, unit):
    rv = np.genfromtxt(os.path.join(dataDir,'misc.dat'), usecols=1)
    return rv*unit


def loadRadius(dataDir, unit, interfaces=False):
    r = np.genfromtxt(os.path.join(dataDir, 'used_rad.dat'))*unit
    dr = r[1:] - r[:-1]
    if not interfaces:
        r = 0.5*(r[1:] + r[:-1])
    return (r, dr)

def loadRadiusFrom1D(dataFile, unit):
    data = np.fromfile(dataFile, dtype=float)
    r = data[::4]*unit
    return r

def loadPhi(dataDir, Nphi=None):
    # Hardcode
    Nphi = Nphi if Nphi is not None else loadNcells(dataDir)[0]
    phi = np.linspace(-np.pi, np.pi, Nphi)
    return phi

def loadMeshGrid(dataDir, unit, interfaces=False):
    # return a meshgrid for the disk to plot data
    Phi, R = loadMeshGridPolar(dataDir, unit, interfaces=interfaces)
    X = R*np.cos(Phi)
    Y = R*np.sin(Phi)
    return (X,Y)

def loadMeshGridPolar(dataDir, unit, interfaces=False):
    # return a meshgrid for the disk to plot data
    phi = loadPhi(dataDir)
    r, dr = loadRadius(dataDir, unit, interfaces=interfaces)
    Phi, R = np.meshgrid(phi, r)
    return (Phi,R)

def loadUnits(dataDir):
    ### load data units
    try:
        units = {l[0] : float(l[1])*u.Unit(l[2]) for l in
                 [l.split() for l in open(os.path.join(dataDir,'units.dat'),'r')
                  if l.split()[0] != '#' and len(l.split())==3]}
        ### fix temperature unit
        units['temperature'] = 1*u.K
    except FileNotFoundError:
        # Fall back to dimensionless units
        units = { bu : 1 for bu in ['mass', 'time', 'length'] }
    return units

def loadNcells(dataDir):
    Nr = len(np.genfromtxt(os.path.join(dataDir, 'used_rad.dat')))-1
    Nphi = int(os.path.getsize(
        os.path.join(dataDir,'gasdens0.dat'))/(8*Nr))
    return (Nphi, Nr)

def load1dRadial(n, dataFilePattern, unit, lengthunit=None):
    data = np.fromfile(dataFilePattern.format(n), dtype=float)
    if 'torque_planet' in dataFilePattern:
        v = data[1::2]*unit
        r = data[::2]
    else:
        v = data[1::4]*unit
        r = data[::4]
    if lengthunit is None:
        return v
    else:
        r = r*lengthunit
        return (r,v)

def load1dRadialAveragedFrom2d(n, dataFilePattern, Nr, Nphi, unit):
    # Load 2d data and average over the azimuthal domain
    data = load2d(n, dataFilePattern, Nr, Nphi, unit)
    rv = np.mean(data, axis=1)
    return rv

def load2d(n, dataFilePattern, Nr, Nphi, unit):
    # Load 2d data an reshape it
    if any([s in dataFilePattern for s in ['vrad']]):
        Nr = Nr+1
    rv = np.fromfile(dataFilePattern.format(n)).reshape(Nr, Nphi)*unit
    return rv

def loadScalar(filepath, varname):
    # check whether the text file is self documented
    found_variables = {}
    with open(os.path.join(filepath)) as f:
        for line in f:
            line = line.strip()
            if line[0] != "#":
                break
            identifier = "#variable:"
            if line[:len(identifier)] == identifier:
                col, name, unitstr = [s.strip() for s in line[len(identifier):].split("|")]
                found_variables[name] = [col, unitstr]
    # get data
    col = found_variables[varname][0]
    unit = u.Unit(found_variables[varname][1])
    data = np.genfromtxt(filepath, usecols=int(col))*unit
    # get time
    timename = "physical time"
    col = found_variables[timename][0]
    unit = u.Unit(found_variables[timename][1])

    time = np.genfromtxt(filepath, usecols=int(col))*unit
    return (time, data)

def loadScalarwithPlanets(dirpath, n, varname):
    # first try to get the variable from the quantities file
    if n == "":
        rv = loadScalar(os.path.join(dirpath, "Quantities.dat"), varname)
        return rv
    else:
    # then try to get it from the planet files if n is set
        rv = loadScalar(os.path.join(dirpath, "bigplanet{}.dat".format(n)), varname)
        return rv

def loadScalarOld(datadir, n, key, varDict, units):
    try:
        # apply the number n to the file path to handle planet files
        datafile = os.path.join( datadir, varDict[key]['file'].format(n))
        datacol = varDict[key]['datacol']
        timecol = varDict[key]['timecol']
        unitdict = varDict[key]['unitpowers']
        t = np.genfromtxt( datafile, usecols=timecol)*units['time']
        v = np.genfromtxt( datafile, usecols=datacol)
        # apply unit to data
        for unit, power in unitdict.items():
            v = v*units[unit]**power
    except KeyError:
        raise KeyError("Don't know how to load variable '{}'".format(key))
    return t,v

class Reader:

    def __init__(self, dataDir):
        self.dataDir = dataDir
        self.units = loadUnits(dataDir)
        # get times of the coarse output steps
        self.outputTimes = np.genfromtxt(os.path.join(dataDir,'misc.dat'))[:,1]*self.units['time']
        self.Nphi, self.Nr = loadNcells(dataDir)
        self.phi = loadPhi(self.dataDir)
        self.r, self.dr = loadRadius(self.dataDir, self.units['length'])
        self.Phi, self.R = loadMeshGridPolar(dataDir, self.units['length'])
        self.X, self.Y = loadMeshGrid(dataDir, self.units['length'])


    def getScalar(self, key, n="", varDict = scalarVars, frame=None, fix_length=True):
        if key in aliases:
            key = aliases[key]

        # first try to load from self documented file
        try:
            t,v = loadScalarwithPlanets(self.dataDir, n, key)
        except KeyError:
            print("scalar {} could not be loaded with new loader".format(key))
            t,v = loadScalarOld(self.dataDir, n, key, varDict, self.units)

        if fix_length and (len(t) != len(v)):
            N = min(len(t), len(v))
            t = u.Quantity(t[:N])
            v = u.Quantity(v[:N])
        # return only the data point that belongs to frame
        if frame is not None:
            tframe = self.outputTimes[frame]
            ind = np.argmin( np.abs(t - tframe))
            return (t[ind],v[ind])
        return (t,v)

    def getScalarPlanet(self, key, n, frame=None):
        offset = scalarVarsPlanet[key]['offset']
        return self.getScalar(key, n=n+offset, varDict=scalarVarsPlanet, frame=frame)

    def get1D(self, key, unit=None):
        # return a TimeSeries1D object which can be indexed via the [] operator
        # and return a list for each timestep containing
        # ['time', 'radius', 'data'] entries
        if unit is None:
            unit = units1D[key]
        ts = TimeSeries1D( self.dataDir, key, unit, self.outputTimes,
                           self.units['length'], name=key)
        return ts

    def load2d(self, n, dataFilePattern, unit):
        return load2d(n, os.path.join(self.dataDir, dataFilePattern), self.Nr, self.Nphi, unit)

class TimeSeries1D:
    # A class to handle on demand loading of 1D data

    def __init__(self, datadir, basename, dataunit, times, lengthunit, name=None):
        self.basename = basename
        self.dataunit = dataunit
        self.lengthunit = lengthunit
        self.times = times
        self.name = basename if name is None else name

        self.datafiles = sorted([os.path.join(datadir,f) for
                                 f in os.listdir(datadir) if basename+'1D' in f],
                           key=lambda x: int(os.path.basename(x)[len(basename)+2:-4]))

        self.N = len(self.datafiles)

    def __getitem__(self, key):
        rv = self._load(key)
        if hasattr(self, '_process'):
            rv = self._process(rv)
        return rv

    def _load(self, key):
        t = self.times[key]
        # load data
        data = np.fromfile(self.datafiles[key], dtype=float)
        r = data[::4]*self.lengthunit
        v = data[3::4]*self.dataunit
        if isinstance(v, u.Quantity):
            v = v.decompose()
        return [t, r, v]
