import numpy as np
import os
import logging
from functools import lru_cache
from astropy import units as u


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
    logging.debug(found_variables)
    return found_variables

def load_text_data_variables(filepath):
    return _load_text_data_variables(filepath, os.path.getmtime(filepath))


@lru_cache(20)
def _load_data(filepath, timestamp):
    return np.genfromtxt(filepath).T

def load_data(filepath):
    return _load_data(filepath, os.path.getmtime(filepath))


def load_text_data_file(filepath, varname, Nmax=np.inf):
    # get data
    variables = load_text_data_variables(filepath)
    col = int(variables[varname][0])
    unit_str = variables[varname][1]
    unit_str = unit_str.replace("1/s", "s-1")
    unit = u.Unit(unit_str)
    file_data = load_data(filepath)
    data = file_data[col] * unit
    if data.isscalar:
        data = u.quantity.Quantity([data])
    return data