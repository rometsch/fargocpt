#!/usr/bin/env python3

import sys
import os
import argparse
import yaml
import copy
import re

old_to_new = {
'ALPHAVISCOSITY': {'newname': 'ViscousAlpha'},
 'Adiabatic': {'newname': 'none', 'hint': 'EquationOfState: Ideal'},
 'AlphaThreshold': {'newname': 'none'},
 'CoolingRadiativeLocal': {'newname': 'none',
  'hint': 'SurfaceCooling: thermal'},
 'CoolingScurve': {'newname': 'none', 'hint': 'SurfaceCooling: scurve'},
 'DT': {'newname': 'MonitorTimestep'},
 'DebugOutputs': {'newname': 'none'},
 'DomegaDrZero': {'newname': 'none', 'hint': 'OuterBoundaryAzi = zeroshear'},
 'ExplicitViscosity': {'newname': 'none',
  'hint': 'Set ViscousAlpha > 0 or ConstantViscosity > 0'},
 'ForcedCircular': {'newname': 'none'},
 'HeatingStar': {'newname': 'none',
  'hint': 'Set temperature of any Nbody object > 0 to enable irradiation from this object.'},
 'HeatingStarFactor': {'newname': 'none',
  'hint': 'Change the temperature of the Nbody object to modulate the heating.'},
 'HeatingStarRampingTime': {'newname': 'none',
  'hint': "Use key 'irradiation ramp-up time' for Nbody objects."},
 'HeatingStarSimple': {'newname': 'none',
  'hint': 'This is now the only option.'},
 'IntegratePlanets': {'newname': 'none'},
 'NSEC': {'newname': 'Naz'},
 'NTOT': {'newname': 'Nsnapshots'},
 'Ninterm': {'newname': 'Nmonitor'},
 'StarRadius': {'newname': 'none',
  'hint': "Use key 'radius' for Nbody objects."},
 'StarTemperature': {'newname': 'none',
  'hint': "Use key 'temperature' for Nbody objects."},
 'StellarRotation': {'newname': 'none',
  'hint': 'InnerBoundaryVaziKeplerianFactor, InnerBoundaryVazi = keplerian'},
 'StsNu': {'newname': 'none', 'hint': 'The STS module has been removed.'},
 'Temperaturecgs0': {'newname': 'none',
  'hint': "Use Temperature0 instead. Append a 'K' to the value to set the temperature in Kelvin."},
 'VISCOSITY': {'newname': 'ConstantViscosity'},
 'VRadIn': {'newname': 'none',
  'hint': 'InnerBoundaryVradKeplerianFactor, InnerBoundaryVrad = keplerian'},
 'ViscosityInCGS': {'newname': 'none',
  'hint': "Use the 'ConstantViscosity' key instead. It supports specifying units."},
 'discmass': {'newname': 'DiskMass'},
 'massoverflow': {'newname': 'RocheLobeOverflow'},
 'mofaveragingtime': {'newname': 'ROFAveragingTime'},
 'mofgamma': {'newname': 'ROFGamma'},
 'mofplanet': {'newname': 'ROFPlanet'},
 'moframpingtime': {'newname': 'ROFRampingTime'},
 'moftemperature': {'newname': 'ROFTemperature'},
 'mofvalue': {'newname': 'ROFValue'},
 'variableTransfer': {'newname': 'ROFVariableTransfer'},
 'zbufferMaxAngle': {'newname': 'none'},
 'zbufferSize': {'newname': 'none'}
}
old_keys = {k.lower(): k for k in old_to_new.keys()}

def replace_word(word, replacement, string):
    pattern = r"(^|\s)" + re.escape(word) + r"\s*:"
    return re.sub(pattern, replacement + ":", string)

def get_new_lines(line, old, new, verbose=False):
    new_lines = []
    if new["newname"] == "none":
        text = "# " + line
        new_lines.append(text)
        if verbose:
            print(text.strip())
        if "hint" in new:
            text = "# hint: " + new["hint"]
            if verbose:
                print(text.strip())                        
            new_lines.append(text + "\n")
        else:
            text = "# has beed removed without replacement"
            new_lines.append(text + "\n")
            if verbose:
                print(text.strip())
    else:
        text = replace_word(old, new["newname"], line)
        new_lines.append(text)
        if verbose:
            print("<"*3)
            print(line.strip())
            print("-"*len(line))
            print(text.strip())
            print(">"*3)
            print()
    return new_lines

def replace_parameter_names(yaml_file, dry=False, verbose=False):
    with open(yaml_file, "r") as infile:
        config = yaml.safe_load(infile)

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

    new_lines = []
    for i, line in enumerate(lines):
        found_line = False
        lline = line.lower()
        for old, new in keys_update.items():
            lold = old.lower()
            if re.search(re.escape(lold) + r"\s*:", lline):
                found_line = True
                new_lines += get_new_lines(line, old, new, verbose=verbose)
                break
        if not found_line:
            new_lines.append(line)

    if not dry:
        with open(yaml_file, "w") as outfile:
            outfile.writelines(new_lines)

def main():
    parser = argparse.ArgumentParser(description='Replace parameter names in a yaml file.')
    parser.add_argument('filename', help='YAML file to process')
    parser.add_argument('--dry', action='store_true', help='Dry run, do not modify file')
    parser.add_argument('--verbose', action='store_true', help='Verbose output. Show replacements.')
    opts = parser.parse_args()

    replace_parameter_names(opts.filename, dry=opts.dry, verbose=opts.verbose)


if __name__ == '__main__':
    main()