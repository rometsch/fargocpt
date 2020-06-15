#!/usr/bin/env python3
# Convert an ini config to a json file.
import os
import json
import argparse


def main():
    args = parse_cli_args()
    params = parse_ini_file(args.infile, comments_enabled=not args.no_comments)
    try:
        if args.no_comments:
            planet_config_path = params["PlanetConfig"]
        else:
            planet_config_path = params["PlanetConfig"]["value"]
        if os.path.exists(planet_config_path):
            planet_params = parse_planet_config(planet_config_path)
            params["Nbody"] = planet_params
        else:
            print("Planet config file '{}' not found.".format(planet_config_path))
            print("Hint: are you in the same dir as the fargo binary?")
    except KeyError:
        pass
    write_json_file(params, args.outfile)


def parse_planet_config(planet_config_file):
    """ Parse a FargoCPT planet config file.

    Parameters
    ----------
    planet_config_file: str
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


def parse_ini_file(file_path, comments_enabled=True):
    params = {}
    comment_counter = 0
    with open(file_path, "r") as in_file:
        for n, line in enumerate(in_file):
            data = parse_line(line, n)
            if len(data) == 0:
                continue
            if data["type"] == "comment":
                if not comments_enabled:
                    continue
                params["comment{}".format(comment_counter)] = data["comment"]
                comment_counter += 1
            elif data["type"] == "value":
                key = data["key"]
                value = data["value"]
                if comments_enabled:
                    del data["type"]
                    del data["key"]
                    params[key] = data
                else:
                    params[key] = value
                # params[key] = value
                # if "comment" in data:
                #     params[key+"_comment"] = data["comment"]
    return params


def write_json_file(params, out_file_path):
    """ Write the parameters to a json file.

    Parameters
    ----------
    params: dict
        Dictionary containing the values.
    out_file_path
        Path to the output file.
    """
    with open(out_file_path, "w") as out_file:
        json.dump(params, out_file, indent=4)


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
    parser.add_argument("outfile", help="Output json file.")
    parser.add_argument("-nc", "--no-comments", default=False,
                        action="store_true", help="Disable all comments.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
