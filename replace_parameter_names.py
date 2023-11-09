#!/usr/bin/env python3

import sys
import os
import argparse
import ruamel.yaml # for handling comments
import copy
import re

old_to_new = {
    "Viscosity" : "ConstantViscosity",
    "AlphaViscosity" : "ViscousAlpha"
}
old_keys = {k.lower(): k for k in old_to_new.keys()}

def replace_word(word, replacement, string):
    pattern = r"(^|\s)" + re.escape(word) + r"\s*:"
    return re.sub(pattern, replacement + ":", string)
    
def replace_parameter_names(yaml_file, dry=False, verbose=False):
    yaml = ruamel.yaml.YAML()
    with open(yaml_file, "r") as infile:
        config = yaml.load(infile)

    # Replace parameter names
    lkeys = {k.lower(): k for k in config.keys()}
    
    keys_update = {}

    for old, new in old_to_new.items():
        lold = old.lower()
        if lold in lkeys:
            keys_update[lkeys[lold]] = new
            
    print("Keys to be updated:", keys_update)
    
    lines = []
    with open(yaml_file, "r") as infile:
        lines = infile.readlines()

    for i, line in enumerate(lines):
        for old, new in keys_update.items():
            if old in line:
                lines[i] = replace_word(old, new, line)
                if verbose:
                    print("<"*3)
                    print(line.strip())
                    print("-"*len(line))
                    print(lines[i].strip())
                    print(">"*3)
                    print()

    if not dry:
        with open(yaml_file, "w") as outfile:
            outfile.writelines(lines)

def main():
    parser = argparse.ArgumentParser(description='Replace parameter names in a yaml file.')
    parser.add_argument('filename', help='YAML file to process')
    parser.add_argument('--dry', action='store_true', help='Dry run, do not modify file')
    parser.add_argument('--verbose', action='store_true', help='Verbose output. Show replacements.')
    opts = parser.parse_args()

    replace_parameter_names(opts.filename, dry=opts.dry, verbose=opts.verbose)


if __name__ == '__main__':
    main()