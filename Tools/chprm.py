#!/usr/bin/env python3
# Replace a parameter value in a text based config file.
# The only assumption is that the parameter name and the value
# are separated by whitespaces and followed by whitespaces or newlines.
import re
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="config file to change parameter in")
    parser.add_argument("param", help="parameter to be changed")
    parser.add_argument("value", help="new value")
    parser.add_argument("-o", help="output file path")
    args = parser.parse_args()

    change_param(args.infile, args.param, args.value,
                 outfile_path = args.o)

def change_param(infile_path, param, value, outfile_path=None):

    outfile_path = outfile_path if outfile_path is not None else infile_path

    ptrn = re.compile(r"\A([ \t]*{}[ \t]+)([\S]+(?=\s))".format(param))

    NreplTot = 0
    with open(infile_path, 'r') as infile:
        outlines = []
        for line in infile:
            out, Nrepl = re.subn(ptrn, r"\g<1>{}".format(value), line)
            outlines.append( out )
            NreplTot += Nrepl

    if NreplTot > 1:
        raise AssertionError("Replaced {} instead of 1 occurances of parameter '{}'".format(Nrepl, param))

    if NreplTot == 0:
        raise AssertionError("Could not find parameter '{}' in '{}'".format(param, infile_path))

    with open(outfile_path, 'w') as outfile:
        outfile.write("".join(outlines))


if __name__=="__main__":
    main()
