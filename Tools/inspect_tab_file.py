#!/usr/bin/env python3
import argparse
import numpy as np
from subprocess import run

def main():
    opt = parse_opt()
    names = parse_header(opt.filename)
    
    if len(opt.column) == 0 or opt.list:
        print_variable_info(names)
    
    if len(opt.column) > 0:
        print_columns(opt.filename, opt.column, formatstrs=opt.format)
    
def print_columns(filename, columns, formatstrs=None):
    data = np.genfromtxt(filename, usecols=columns)
    
    autoformat = False
    if formatstrs is None:
        formatstrs = [".3g"] * len(columns)
        autoformat = True
    s = ""
    for dataline in data:
        colums = []
        for fs, x in zip(formatstrs, dataline):
            if autoformat and x.is_integer():
                fs = "{}"
                x = int(x)
            else:
                fs = "{:" + fs + "}"
            colums.append("{:10s}".format(fs.format(x)))
        s += "\t".join(colums) + "\n"

    run("less", text=True, input=s)

    
def parse_header(filename):
    header = get_header(filename)
    names = parse_syntax_line(header)
    return names

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
    parser.add_argument("column", type=int, nargs="*", default=[], help="Column(s) to inspect.")
    parser.add_argument("--list", action="store_true", help="List the available variables if possible.")
    parser.add_argument("--format", type=str, nargs="*", default=None, help="Format strings for the columns.")
    opt = parser.parse_args()
    return opt


if __name__=="__main__":
    main()