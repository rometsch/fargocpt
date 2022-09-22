#!/usr/bin/env python3
import argparse
import os
import numpy as np
from subprocess import run


def main():
    opt = parse_opt()
    names, units = parse_header(opt.filename)
    if opt.units is None:
        units = {}
        tounits = []
    else:
        if len(opt.units) != 0 and len(opt.column) != len(opt.units):
            raise ValueError("Arguments of units must be the same number as columns!")
        tounits = opt.units

    if len(opt.column) == 0 or opt.list:
        print_variable_info(names)

    if len(opt.column) > 0:
        data = load_data(opt.filename, opt.column, names, units, tounits)
        print_columns(data, formatstrs=opt.format)


def load_data(filename, columns, names, units, tounits):
    if len(units) > 0:
        import astropy.units as u
    data = np.genfromtxt(filename, usecols=columns, unpack=True)
    rv = []
    for n, k in enumerate(columns):
        name = names[k]
        x = data[n]
        if name in units:
            x = x*u.Unit(units[name])
            x = x.decompose()
            if len(tounits) > 0:
                x = x.to(tounits[n])
                print(x)
        rv.append(x)
    return rv

def print_columns(data, formatstrs=None):
    Nvars = len(data)
    autoformat = False
    if formatstrs is None:
        formatstrs = [".3g"] * Nvars
        autoformat = True
    output_str = ""
    for n in range(len(data[0])):
        values = [x[n] for x in data]
        colums = []
        for fs, x in zip(formatstrs, values):
            unit = ""
            if hasattr(x, "unit"):
                unit += f" {x.unit}"
                val = x.value
            else:
                val = x
            if autoformat and val.is_integer():
                fs = "{}"
                val = int(val)
            else:
                fs = "{:" + fs + "}"
            s = "{:>10s}".format(fs.format(val)) + unit
            colums.append(s)
            
        output_str += "\t".join(colums) + "\n"

    run("less", text=True, input=output_str)


def parse_header(filename):
    header = get_header(filename)
    detailed_header = False
    for line in header:
        line = line.strip()
        if line.lower().startswith("variable:"):
            detailed_header = True

    units = {}
    if detailed_header:
        names, units = parse_detailed_header(header)
    else:
        names = parse_syntax_line(header)
    return names, units


def parse_detailed_header(lines):
    names = []
    units = {}
    for line in lines:
        if line.lower().startswith("variable:"):
            name = line.split(":", 1)[1]
            unit = name.split("|")[2].strip()
            name = name.split("|")[1].strip()
            names.append(name)
            units[name] = unit
    return names, units


def print_variable_info(names):
    print("Available variables:")
    for n, name in enumerate(names):
        print("{:2}   {}".format(n, name))


def get_header(filename):
    lines = []
    with open(filename, "r") as infile:
        for line in infile:
            line = line.strip()
            if line[0] == "#":
                pass
            elif line == "":
                continue
            else:
                break
            line = line.lstrip("#").lstrip()
            lines.append(line)
    return lines


def parse_syntax_line(header_lines):
    syntax_line = None
    for line in header_lines:
        if line.lower().startswith("syntax"):
            syntax_line = line
    if syntax_line is not None:
        names_string = syntax_line.split(":", 1)[1]
    else:
        return None
    names = [n.strip() for n in names_string.split("<tab>")]
    return names


def parse_opt():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="Path of the file to be inspected.")
    parser.add_argument("column", type=int, nargs="*",
                        default=[], help="Column(s) to inspect.")
    parser.add_argument("--list", action="store_true",
                        help="List the available variables if possible.")
    parser.add_argument("--format", type=str, nargs="*",
                        default=None, help="Format strings for the columns.")
    parser.add_argument("-u", "--units", nargs="*", type=str, default=None,
                        help="Attempt to load units and apply them.")
    opt = parser.parse_args()
    return opt


if __name__ == "__main__":
    main()
