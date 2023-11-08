import numpy as np
import astropy.units as u
import re
import yaml

from read_par_file import read_unit_file

class read1D:
    def __init__(self, output_folder_path, quantity):
        self.output_folder_path = output_folder_path
        self.quantity = quantity

        self.units = read_unit_file(output_folder_path)
        self.l0_to_cm = self.units["length"]

        with open(self.output_folder_path + "/info1D.yml") as f:
            data = yaml.safe_load(f)
            info = data[quantity]            
            self.Nr = int(info["Nrad"])
            self.unit = u.Unit(info["unit"])
            # self.code_to_cgs_factor = float(info["code_to_cgs_factor"])
            # print(quantity, self.unit, self.code_to_cgs_factor)

    def read(self, dt, return_min_max=False):
        with open(self.output_folder_path + f"/snapshots/{dt}/" + str(self.quantity) + "1D.dat", "rb") as f:
            Nr = self.Nr
            data = np.fromfile(f)
            rs = data[::4]
            avg_values = data[1::4]

            if return_min_max:
                min_value = data[2::4]
                max_value = data[3::4]

            if len(avg_values) != Nr:
                raise ValueError("read1D.py: Output data length ({}) is not a multiple of Nr = {}".format(len(avg_values), Nr))

            avg_values = avg_values*u.Unit(self.unit)
            rs = rs*u.Unit(self.units["length"])

            if return_min_max:
                min_value *= u.Unit(self.unit)
                max_value *= u.Unit(self.unit)


            if return_min_max:
                return rs, avg_values, min_value, max_value
            else:
                return rs, avg_values


if __name__ == "__main__":
    mass_dens_reader = read1D("/home/jordan/develop/fargocpt/test/gasMassFlow1D/out/", 'Sigma')

    rs, density = mass_dens_reader.read(0)
    print(rs)
    print(density)
