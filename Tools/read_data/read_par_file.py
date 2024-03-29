import sys
import os

import numpy as np
import astropy.units as u

import yaml

def read_par_file(path_to_file):
    """
    input:  path the the *.par file
    return: dictionary with the content of the par file
    """
    unit_dict = {}
    if os.path.isfile(path_to_file):
        with open(path_to_file) as f:
            units = f.readlines()

        for line in units:
            l = line.replace('\t', ' ')
            l = l.replace('\n', '')
            l = l.split('#', 1)[0]
            l = l.split(' ')
            l = list(filter(None, l))
            if len(l) > 1:
                if l[1][0].isdigit():
                    if l[1].isdigit():
                        unit_dict[l[0]] = int(l[1])
                    else:
                        tmp = l[1].split(',')
                        if len(tmp) == 1:
                            unit_dict[l[0]] = float(l[1])
                        else:
                            unit_dict[l[0]] = tmp
                elif l[1].lower() == 'yes':
                        unit_dict[l[0]] = True
                elif l[1].lower() == 'no':
                        unit_dict[l[0]] = False
                else:
                        unit_dict[l[0]] = l[1]
    else:
        raise ValueError("Error in read_par_file.py, could not find : " + path_to_file + "!\n")

    unit_dict['R'] = 1
    unit_dict['gamma'] = unit_dict['AdiabaticIndex']
    unit_dict['alpha'] = unit_dict['AlphaViscosity']

    return unit_dict




def read_unit_file(path_to_folder):

    units_file_name = path_to_folder + "/units.yml"

    unit_dict = {}
    if os.path.isfile(units_file_name):
        with open(units_file_name) as f:
            data = yaml.safe_load(f)
        unit_dict = {}
        for key, val in data.items():
            unit_dict[key] = u.Unit(val["unit"])
    else:
        raise ValueError("Error in read_units, could not find : " + units_file_name  + "!\n")

    return unit_dict



if __name__ == "__main__":
    read_unit_file("/home/jordan/UnitTest/fargocpt/Tools/test_problems/gasMassFlow1D/out")
