import numpy as np
import astropy.units as u
import re
import matplotlib.pyplot as plt


class gasMassFlow1D:
    def __init__(self, output_folder_path, print_unit="solMass/yr"):
        self.output_folder_path = output_folder_path
        with open(self.output_folder_path + "MassFlow1D.info") as f:
            header = f.readline()
            header = f.readline()
            header = f.readline()
            header = f.readline()

            self.Nr = int(re.search("Nr = ([\d]*)", header).groups()[0])

            header = f.readline()
            self.unit = re.search("unit = ([^,]*)", header).groups()[0].strip()
            self.print_unit = print_unit

        self.DT = np.loadtxt(self.output_folder_path + "/monitor/Quantities.dat", usecols=1, dtype=np.int)
    def read(self, dt):
        with open(self.output_folder_path + f"/snapshots/{dt}/MassFlow1D.dat", "rb") as f:
            Nr = self.Nr
            data = np.fromfile(f).reshape(Nr, 4)
            rs = data[:,0]
            data = data[:,1]
            data = data*u.Unit(self.unit).to(self.print_unit)
            return rs, self.DT, data


if __name__ == "__main__":
    mass_flow_reader = gasMassFlow1D("/home/jordan/develop/fargocpt/test/gasMassFlow1D/out/")

    rs, DT, data = mass_flow_reader.read(0)
    print(data)
