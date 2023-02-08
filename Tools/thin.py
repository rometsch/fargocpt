#!/usr/bin/env python3
"""Thin out the output of a fargocpt simulation. E.g. keep only every 100th snapshot.
"""

import os
import argparse
import numpy as np
from pathlib import Path


def main():
    opts = parse_args()
    
    to_inspect = not hasattr(opts, "N")
    if to_inspect:
        inspect(opts.outputdir)
    else:
        thin(opts.outputdir, opts.N, force=opts.y)

def thin(outputdir, N, force=False):
    inds = get_inds(outputdir)
    print(f"Output at path {outputdir}")
    new_inds = inds[::N]
    new_inds = np.append(new_inds, inds[-1])
    print(f"Thinning down to {len(new_inds)} snapshots.")
    print(new_inds)
    print_size(outputdir, inds[0], len(new_inds))

def inspect(outputdir):
    """Print out how many snapshots the directory contains.

    Args:
        outputdir (str): Path to the output directory.
    """
    inds = get_inds(outputdir)
    print(f"Output at path {outputdir}")
    print(f"contains {len(inds)} snapshots.")
    
    print_size(outputdir, inds[0], len(inds))

def print_size(outputdir, Nref, Nsnapshots):

    first_snapshot_dir = Path(outputdir) / "snapshots" / f"{Nref}"

    size_snapshot = size_of_dir(first_snapshot_dir)
    size_rest = 0
    for p in Path(outputdir).glob("*"):
        if p.name == "snapshots":
            continue
        size_rest += size_of_dir(p)
    

    print(f"Size of one snapshot: {sizeof_fmt(size_snapshot)}")
    print(f"Size of all snapshot: {sizeof_fmt(Nsnapshots*size_snapshot)}")
    print(f"Size of rest: {sizeof_fmt(size_rest)}")
    print(f"Size total: {sizeof_fmt(size_rest+Nsnapshots*size_snapshot)}")

def size_of_dir(path):
    return sum(f.stat().st_size for f in Path(path).glob('**/*') if f.is_file())
    

def get_inds(outputdir):
    """Look up the inds of snapshots in the directory.

    Args:
        outputdir (str): Path to the output directory.
        
    Returns:
        list of int: List of indices of snapshots.
    """
    list_file = os.path.join(outputdir, "snapshots", "list.txt")
    inds = np.genfromtxt(list_file, dtype=int)
    return inds

def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    inspect_parser = subparsers.add_parser("inspect")
    inspect_parser.add_argument("outputdir", type=str)
    
    thin_parser = subparsers.add_parser("thin")
    thin_parser.add_argument("outputdir", type=str)
    thin_parser.add_argument("N", type=int, help="Keep every Nth snapshot.")
    thin_parser.add_argument("-y", action="store_true", help="Answer yes to all questions.")

    opts = parser.parse_args()
    return opts


def sizeof_fmt(num, suffix="B"):
    """Convert a number of bytes to human readable string.
    
    https://stackoverflow.com/a/1094933
    CC BY-SA 4.0

    Args:
        num (int): Size in bytes.
        suffix (str, optional): unit suffiux

    Returns:
        str: Human readable size.
    """
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return f"{num:3.1f} {unit}{suffix}"
        num /= 1024.0
    return f"{num:.1f} Yi{suffix}"

if __name__=="__main__":
    main()