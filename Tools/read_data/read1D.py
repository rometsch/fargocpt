import numpy as np
import astropy.units as u
import re


class read1D:
    def __init__(self, output_folder_path, quantity):
        self.output_folder_path = output_folder_path + quantity + "1D"
        with open(self.output_folder_path + ".info") as f:
            self.header = f.readline()
            self.header = f.readline()
            self.header = f.readline()
            self.header = f.readline()

            self.Nr = int(re.search("Nr = ([\d]*)", self.header).groups()[0])
            self.header = f.readline()
            self.unit = re.search("unit = ([^,]*)", self.header).groups()[0].strip()

    def read(self, dt, return_min_max=False):
        with open(self.output_folder_path + str(dt) +  ".dat", "rb") as f:
            Nr = self.Nr
            data = np.fromfile(f)
            rs = data[::4]
            density = data[1::4]

            if return_min_max:
                min_density = data[2::4]
                max_density = data[3::4]

            if len(density) != Nr:
                raise ValueError("read1D.py: Output data length ({}) is not a multiple of Nr = {}".format(len(density), Nr))

            density = density*u.Unit(self.unit)
            rs = rs*u.Unit('au')

            if return_min_max:
                min_density *= u.Unit(self.unit)
                max_density *= u.Unit(self.unit)


            if return_min_max:
                return rs, density, min_density, max_density
            else:
                return rs, density


if __name__ == "__main__":
    mass_dens_reader = read1D("/home/jordan/UnitTest/fargocpt/Tools/test_problems/gasMassFlow1D/out/", 'gasdens')

    rs, density = mass_dens_reader.read(0)
