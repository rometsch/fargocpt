#!/usr/bin/env python3
# Convert an ini config to a json file.
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
            planet_params = parse_planet_config(planet_config_path)
            params["planets"] = planet_params
        else:
            print("Planet config file '{}' not found.".format(planet_config_path))
            print("Hint: are you in the same dir as the fargo binary?")
    except KeyError:
        pass
    write_yaml_file(params, args.outfile)

    insert_comments(comments, args.outfile)

def insert_comments(comments, filename):

    with open(filename, "r") as in_file:
        lines = in_file.readlines()

    inline_comments = {}
    standalone_comments = {}
    first_comment = ""

    new_lines = []

    for com in comments:
        if "this_key" in com:
            inline_comments[com["this_key"]] = com["comment"]
        elif "last_key" in com:
            standalone_comments[com["last_key"]] = com["comment"]
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
            if key in standalone_comments:
                new_lines.append("")
                new_lines.append("# " + standalone_comments[key])
                new_lines.append("")

        new_lines.append(new_line)

    with open(filename, "w") as out_file:
        for line in new_lines:
            print(line, file=out_file)

def parse_planet_config(planet_config_file):
    """ Parse a FargoCPT planet config file.

    Parameters
    ----------first_comment
        Path to the planet config file.
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
    planets = []
    with open(planet_config_file, "r") as in_file:
        for line in in_file:
            planet = dict()
            line = line.strip()
            if line == "" or line[0] == "#":
                continue
            values = line.split()
            for key, val in zip(keys, values):
                planet[key] = val
        try:
            del planet["Nbody interaction"]
        except KeyError:
            pass
        planets.append(planet)
    return planets


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
                comment = {"comment" : data["comment"]}
                if last_key != "":
                    comment["last_key"] = last_key
                comments.append(comment)
                
            elif data["type"] == "value":
                key = data["key"]
                value = data["value"]

                params[key] = value
                if "comment" in data:
                    comment = { "comment" : data["comment"],
                                "this_key" : key}
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
        yaml.dump(params, out_file, width=200, sort_keys=None)


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
    parser.add_argument("-c", "--comments", help="Include comments", action="store_true")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()