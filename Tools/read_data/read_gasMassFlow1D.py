import numpy as np
import astropy.units as u
import re
import matplotlib.pyplot as plt


class gasMassFlow1D:
    def __init__(self, output_folder_path, print_unit="solMass/yr"):
        self.output_folder_path = output_folder_path
        with open(self.output_folder_path + "gasMassFlow1D.info") as f:
            header = f.readline()
            header = f.readline()

            self.Nr = int(re.search("Nr = ([\d]*)", header).groups()[0])

            header = f.readline()
            self.unit = re.search("unit = ([^,]*)", header).groups()[0].strip()
            self.print_unit = print_unit

        self.DT = np.loadtxt(self.output_folder_path + "Quantities.dat", usecols=1, dtype=np.int)
    def read(self):
        with open(self.output_folder_path + "gasMassFlow1D.dat", "rb") as f:
            Nr = self.Nr
            data = np.fromfile(f)
            rs = data[:Nr]
            data = data[Nr:]
            data = data
            num_outputs = len(data)/Nr
            if num_outputs%1 != 0:
                raise ValueError("Output data length ({}) is not a multiple of Nr = {}".format(len(data), Nr))
            num_outputs = int(num_outputs)
            data = data.reshape( int(num_outputs), Nr)
            data = data*u.Unit(self.unit).to(self.print_unit)
            return rs, self.DT, data, num_outputs


if __name__ == "__main__":
    mass_flow_reader = gasMassFlow1D("/home/jordan/UnitTest/fargocpt/Tools/test_problems/gasMassFlow1D/out/")

    rs, DT, data, Noutputs = mass_flow_reader.read()
    print(data[0])
