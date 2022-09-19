#!/usr/bin/env python3
# Convert an ini config to a json file.
from email.policy import default
import os
from typing import OrderedDict
import yaml
import argparse


def main():
    args = parse_cli_args()
    params, comments = parse_ini_file(args.infile)
    try:
        planet_config_path = params["PlanetConfig"]
        if os.path.exists(planet_config_path):
            nbody_params = parse_planet_config(planet_config_path)
            comments.append({"comment": "Nbody", "next_key": "nbody"})
            params["nbody"] = nbody_params
            print(f"Added {len(nbody_params)} planets")
            remove_entry(params, "PlanetConfig")
        else:
            print("Planet config file '{}' not found.".format(planet_config_path))
            print("Hint: are you in the same dir as the fargo binary?")
    except KeyError:
        pass

    handle_default_star(params)
    handle_implicit_units(params)
    for planet in params["nbody"]:
        add_unit(planet, "temperature", "K")
        add_unit(planet, "radius", "solRadius")

    remove_deprecated_entries(params)

    write_yaml_file(params, args.outfile)

    insert_comments(comments, args.outfile)


def remove_entry(params, key):
    if contains(params, key):
        params.pop(keyname(params, key))
    print(f"Removed deprecated parameter {key} which has no effect anymore.")


def remove_deprecated_entries(params):
    """ Remove deprecated settings which have no effect anymore.

    Parameters
    ----------
    params: dict
        Dictionary containing the config.
    """
    obsolete = [
        "Sigma0InCodeUnits",
        "ViscosityInCgs",
        "TemperatureCGS0",
        "HeatingStar",
        "HeatingStarFactor",
        "HeatingStarSimple",
        "HeatingStarRampingTime"
    ]
    for key in obsolete:
        remove_entry(params, key)

    for planet in params["nbody"]:
        try:
            if not get_flag(planet, "irradiate"):
                planet["temperature"] = "0 K"
        except KeyError:
            pass

        remove_entry(planet, "Feels Disk")
        remove_entry(planet, "Nbody interaction")

def get_flag(params, key):
    """ Get a flag from the params dict and return true or false.

    Parameter
    ---------
    params: dict
        Dictionary containing the config.
    search_key: str
        Key to search for.

    Return
    ------
    bool
    """
    if contains(params, key):
        val = params[keyname(params, key)]
        if val.lower() in ["no", "false"]:
            return False
        elif val.lower() in ["yes", "true"]:
            return True
        else:
            raise ValueError(f"{val} can't be interpreted as bool")
    else:
        raise KeyError(f"Key '{key}' not found.")

def contains(d, search_key):
    """ Is there a key in d which in lowercase matches the search key.

    Parameter
    ---------
    params: dict
        Dictionary containing the config.
    search_key: str
        Key to search for.

    Return
    ------
    bool
    """
    found_key = False
    for key in [k for k in d.keys()]:
        if key.lower() == search_key.lower():
            # check if setting is set to yes
            found_key = True
    return found_key


def keyname(d, search_key):
    """ Which key in d matches the search key in lowercase form.

    Parameter
    ---------
    params: dict
        Dictionary containing the config.
    search_key: str
        Key to search for.

    Return
    ------
    str or None
    """
    keyname = None
    for key in [k for k in d.keys()]:
        if key.lower() == search_key.lower():
            # check if setting is set to yes
            keyname = key
    return keyname


def add_unit(params, key, unit):
    if contains(params, key):
        if not unit in params[keyname(params, key)]:
            params[keyname(params, key)] += " " + unit
            print(f"Added unit {unit} to parameter {key}.")


def handle_implicit_units(params):
    """ Add units to parameters where they were implicitly assumed before.

    E.g. the surface density was given in g/cm2 before.
    Now this entry would be in code units!

    Parameter
    ---------
    params: dict
        Dictionary containing the config.
    """
    add_unit(params, "sigma0", "g/cm2")
    add_unit(params, "ParticleRadius", "cm")
    add_unit(params, "ParticleDensity", "g/cm3")
    add_unit(params, "MaximumTemperature", "K")
    add_unit(params, "MinimumTemperature", "K")
    add_unit(params, "mofvalue", "solMass/yr")


def handle_default_star(params):
    """ Add default star to planet if needed. 

    The default star setting is deprecated in the new yml config format.
    Thus, add the default star explicitly on conversion.

    Parameter
    ---------
    params: dict
        Dictionary containing the config.
    """
    found_key = False
    for key in [k for k in params.keys()]:
        if key.lower() == "defaultstar":
            # check if setting is set to yes
            found_key = params[key][0].lower() == "y"
            params.pop(key)

    if not found_key:
        return

    if "StarTemperature" in params:
        temperature = params.pop("StarTemperature") + " K"
    else:
        temperature = "5778 K"

    try:
        if not get_flag(params, "HeatingStar"):
            temperature = "0"
    except KeyError:
        pass

    if "StarRadius" in params:
        radius = params.pop("StarRadius") + " solRadius"
    else:
        radius = "1 solRadius"

    default_star = {
        "name": "DefaultStar",
        "semi-major axis": "0.0 au",
        "mass": "1.0",
        "eccentricity": "0.0",
        "radius": radius,
        "temperature": temperature
    }
    try:
        params["nbody"] = [default_star] + params["nbody"]
        print("Added a default star to the nbody list.")
    except Exception: # If no planet file exists, nbody does not exist and we get error
        params["nbody"] = [default_star]
        print("Added a default star to the nbody list.")


def insert_comments(comments, filename):

    with open(filename, "r") as in_file:
        lines = in_file.readlines()

    inline_comments = {}
    standalone_comments_after = {}
    standalone_comments_before = {}

    first_comment = ""

    new_lines = []

    for com in comments:
        if "this_key" in com:
            inline_comments[com["this_key"]] = com["comment"]
        elif "last_key" in com:
            standalone_comments_after[com["last_key"]] = com["comment"]
        elif "next_key" in com:
            standalone_comments_before[com["next_key"]] = com["comment"]
        else:
            first_comment = com["comment"]

    if first_comment != "":
        new_lines.append("# " + first_comment)
        new_lines.append("")

    for line in lines:
        try:
            res = yaml.load(line, Loader=yaml.Loader)
        except Exception:
            print(line)
            continue

        new_line = line.rstrip()

        if isinstance(res, dict):
            key = list(res.keys())[0]
            if key in inline_comments:
                new_line += "   # " + inline_comments[key]
            if key in standalone_comments_before:
                new_lines.append("")
                new_lines.append("# " + standalone_comments_before[key])
                new_lines.append("")

        new_lines.append(new_line)

        if isinstance(res, dict):
            key = list(res.keys())[0]
            if key in standalone_comments_after:
                new_lines.append("")
                new_lines.append("# " + standalone_comments_after[key])
                new_lines.append("")

    with open(filename, "w") as out_file:
        for line in new_lines:
            print(line, file=out_file)


def parse_planet_config(nbody_config_file):
    """ Parse a FargoCPT planet config file.

    Parameters
    ----------first_comment
        Path to the nbody config file.
    """
    keys = [
        "name",
        "semi-major axis",
        "mass",
        "accretion efficiency",
        "feels disk",
        "Nbody interaction",
        "eccentricity",
        "radius",
        "temperature",
        "irradiate",
        "phi",
        "ramp-up time"
    ]
    nbody = []
    with open(nbody_config_file, "r") as in_file:
        for line in in_file:
            planet = dict()
            line = line.strip()
            if line == "" or line[0] == "#":
                continue
            values = line.split()
            for key, val in zip(keys, values):
                planet[key] = val
        nbody.append(planet)
    return nbody


def parse_ini_file(file_path):
    params = dict()
    comments = []
    last_key = ""
    with open(file_path, "r") as in_file:
        lines = in_file.readlines()
        for n, line in enumerate(lines):
            data = parse_line(line, n)
            if len(data) == 0:
                continue
            if data["type"] == "comment":
                comment = {"comment": data["comment"]}
                if last_key != "":
                    comment["last_key"] = last_key
                comments.append(comment)

            elif data["type"] == "value":
                key = data["key"]
                value = data["value"]

                params[key] = value
                if "comment" in data:
                    comment = {"comment": data["comment"],
                               "this_key": key}
                    comments.append(comment)
                last_key = key
    return params, comments


def write_yaml_file(params, out_file_path):
    """ Write the parameters to a yaml file.

    Parameters
    ----------
    params: dict
        Dictionary containing the values.
    out_file_path
        Path to the output file.
    """
    with open(out_file_path, "w") as out_file:
        yaml.dump(params, out_file, width=200, sort_keys=None, default_style='')


def parse_line(line, n):
    """ Parse a single ini file line.

    Parameters
    ----------
    line: str
        Line from an ini file.
    n: int
        Line number for error message.

    Returns
    -------
    dict
        Dict containing the data.
    """
    line = line.strip()
    if line == "":
        return {}
    if line[0] == "#":
        ent = {
            "type": "comment",
            "comment": line.strip().lstrip("#").strip()
        }
        return ent
    parts = line.strip().split(maxsplit=2)
    if len(parts) == 1:
        print("Error: Line {} only has a key and no value:".format(n+1))
        print(line)
        print("Please investigate manually! Bye bye.")
        exit()
    ent = {
        "type": "value",
        "key": parts[0],
        "value": parts[1]
    }
    # handle empty string ("") as value
    if ent["value"] == "\"\"":
        ent["value"] = ""
    if len(parts) == 3:
        ent["comment"] = parts[2].lstrip("#").strip()
    return ent


def parse_cli_args():
    """ Parse command line arguments.

    Return
    ------
    Namespace
        Arguments inside a namespace.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Path of INI file to be parsed.")
    parser.add_argument("outfile", help="Output yaml file.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
