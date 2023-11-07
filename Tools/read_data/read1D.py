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
            self.code_to_cgs_factor = float(info["code_to_cgs_factor"])
            print(quantity, self.unit, self.code_to_cgs_factor)

    def read(self, dt, return_min_max=False):
        with open(self.output_folder_path + f"/snapshots/{dt}/" + str(self.quantity) + "1D.dat", "rb") as f:
            Nr = self.Nr
            data = np.fromfile(f)
            rs = data[::4]
            density = data[1::4]

            if return_min_max:
                min_density = data[2::4]
                max_density = data[3::4]

            if len(density) != Nr:
                raise ValueError("read1D.py: Output data length ({}) is not a multiple of Nr = {}".format(len(density), Nr))

            density = density*u.Unit(self.unit)*self.code_to_cgs_factor
            rs = rs*self.l0_to_cm

            if return_min_max:
                min_density *= u.Unit(self.unit)*self.code_to_cgs_factor
                max_density *= u.Unit(self.unit)*self.code_to_cgs_factor


            if return_min_max:
                return rs, density, min_density, max_density
            else:
                return rs, density


if __name__ == "__main__":
    mass_dens_reader = read1D("/home/jordan/develop/fargocpt/test/gasMassFlow1D/out/", 'Sigma')

    rs, density = mass_dens_reader.read(0)
    print(rs)
    print(density)
