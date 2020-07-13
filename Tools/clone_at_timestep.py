#!/usr/bin/env python3
""" Extract output data for one timestep from a fargocpt simulation begin a new one at that time.
"""
import os
import re
import argparse
from shutil import copy2, copytree, ignore_patterns


def main():
    args = parse_cmd_args()

    create_new_simdir(args.simdir, args.dstdir)
    src_output = os.path.join(args.simdir, args.out_dir_name)
    dst_output = os.path.join(args.dstdir, args.out_dir_name)
    os.makedirs(dst_output)
    copy_binary_data(src_output, dst_output, args.timestep)
    copy_scalar_data(src_output, dst_output, args.timestep)


def create_new_simdir(src_dir, dst_dir, ignore=None):
    """ Create a new simulation dir to copy data to.

    Parameters
    ----------
    src_dir: str
        Path to the source simulation directory.
    dst_dir: str
        Path to the destination simulation directory.
    ignore: list of str
        Names to ignore when copying the dir.
    """
    if ignore is None:
        ignore = ["output", "outputs", "plots",
                  "core.*", "out.log*", "err.log"]
    copytree(src_dir, dst_dir, ignore=ignore_patterns(*ignore))


def copy_binary_data(src_dir, dst_dir, timestep):
    """ Copy binary output data to the new directory. 

    Parameters
    ----------
    src_dir: str
        Path to the source output directory.
    dst_dir: str
        Path to the destination output directory.
    timestep: int
        Timestep to be extracted.
    """
    variables = get_binary_variable_list(src_dir)
    for var in variables:
        src_filename = "{}{}.dat".format(var, timestep)
        src = os.path.join(src_dir, src_filename)
        dst_filename = "{}{}.dat".format(var, 0)
        dst = os.path.join(dst_dir, dst_filename)
        copy2(src, dst)
        try:
            copy2(os.path.join(src_dir, f"{var}.info"), dst_dir)
        except FileNotFoundError:
            pass


def copy_scalar_data(src_dir, dst_dir, timestep):
    """ Copy scalar output data to the new directory. 

    Parameters
    ----------
    src_dir: str
        Path to the source output directory.
    dst_dir: str
        Path to the destination output directory.
    timestep: int
        Timestep to be extracted.
    """

    copy_scalar_data_single_timestep(
        os.path.join(src_dir, "Quantities.dat"),
        os.path.join(dst_dir, "Quantities.dat"),
        timestep)
    copy_scalar_data_single_timestep(
        os.path.join(src_dir, "misc.dat"),
        os.path.join(dst_dir, "misc.dat"),
        timestep)
    for n in range(1000):
        try:
            file_name = f"planet{n}.dat"
            copy_scalar_data_single_timestep(
                os.path.join(src_dir, file_name),
                os.path.join(dst_dir, file_name),
                timestep)
            file_name = f"bigplanet{n}.dat"
            copy_scalar_data_single_timestep(
                os.path.join(src_dir, file_name),
                os.path.join(dst_dir, file_name),
                timestep)
        except FileNotFoundError:
            pass


def copy_scalar_data_single_timestep(src_file, dst_file, timestep):
    """ Copy a scalar text based output file only keeping a single data line.

    Parameters
    ----------
    src_file: str
        Path to the source data file.
    dst_file: str
        Path to the destination data file.
    timestep: int
        Timestep for which the line is to be extracted.
    """
    lines = extract_lines_scalar_output(src_file, timestep)
    with open(dst_file, "w") as out_file:
        for line in lines:
            print(line, file=out_file)


def extract_lines_scalar_output(file_path, timestep):
    """ Extract a single data line from a text based output file.

    Use the first occurance that matches the requirement.
    Copy the header.

    Parameters
    ----------
    file_path: str
        Path to the data file.
    timestep: int
        Timestep for which the line is to be extracted.
    """
    out_lines = []
    with open(file_path, "r") as in_file:
        for line in in_file:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == "#":
                out_lines.append(line)
                continue
            if line.split()[0] == f"{timestep}":
                line = zero_step_counters(line)
                out_lines.append(line)
                break
    return out_lines

def zero_step_counters(line):
    """ Replace step counters by 0.
    
    Fargocpt outputs these step counters in the first and sometimes second column.
    
    Parameters
    ----------
    line: str
        Line in which to replace step counters by 0.
        
    Returns
    -------
    str
        Modified line.
    """
    second = line.split()[1]
    if re.match(r"^\d*$", second):
        return "\t".join(["0", "0"] + line.split()[2:])
    else:
        return "\t".join(["0"] + line.split()[1:])
    

def get_binary_variable_list(output_dir):
    """ Get a list of variables that were outputted.

    When restarting, fargocpt determines which 1D files to write depending
    on what is already writter out.

    Parameters
    ----------
    output_dir: str
        Path to the output directory.

    Returns
    -------
    list of str
        A list of variable names.
    """
    files = os.listdir(output_dir)
    ptrn = re.compile(r"(gas.*[\D])[0-9]+\.dat")
    names = set()
    for f in files:
        res = re.match(ptrn, f)
        if res is not None:
            names.add(res.groups()[0])
    # fix inconsitent naming of gastau2 var
    if "gastau" in names:
        names.remove("gastau")
        names.add("gastau2")
    return names


def parse_cmd_args():
    """ Parse command line arguments.

    Returns
    -------
    namespace
        Namespace containing the key value pairs.
    """
    desc = "Extract data for a single timestep from a fargocpt simulation."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("simdir", help="Source simulation dir.")
    parser.add_argument("dstdir", help="Path of the new simulation dir.")
    parser.add_argument("timestep", type=int, help="Timestep to be extracted.")
    parser.add_argument("--out-dir-name", default="output",
                        help="Name of the output dir inside the simulation dir.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
