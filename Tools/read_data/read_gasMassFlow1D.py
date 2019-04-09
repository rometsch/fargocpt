import numpy as np
import astropy.units as u
import re
import matplotlib.pyplot as plt

with open("gasMassFlow1D.info") as f:
    header = f.readline()
    Nr = int(re.search("Nr = ([\d]*)", header).groups()[0])
    unit = re.search("unit = ([^,]*)", header).groups()[0].strip()

with open("gasMassFlow1D.dat", "rb") as f:
    data = np.fromfile(f)
    rs = data[:Nr]
    data = data[Nr:]
    Noutputs = len(data)/Nr
    if Noutputs%1 != 0:
        raise ValueError("Output data length ({}) is not a multiple of Nr = {}".format(len(data), Nr))
    Noutputs = int(Noutputs)
    data = data.reshape( int(Noutputs), Nr)
    data = data*u.Unit(unit).to("solMass/yr")

fig, ax = plt.subplots()             
for n in range(Noutputs):
    ax.plot(rs, data[n], label="{}".format(n))
plt.legend()                           
plt.show()
