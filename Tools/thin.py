#!/usr/bin/env python3
"""Thin out the output of a fargocpt simulation. E.g. keep only every 100th snapshot.
"""

import os
import argparse
import numpy as np
from pathlib import Path
import shutil


def main():
    opts = parse_args()
    
    to_inspect = not hasattr(opts, "N")
    if to_inspect:
        inspect(opts.outputdir)
    else:
        thin(opts.outputdir, opts.N, force=opts.y)

def thin(outputdir, N, force=False):
    outputdir = Path(outputdir)
    inds = get_inds(outputdir)
    print(f"Output at path {outputdir}")
    
    tmpfile = outputdir/"thin_newinds.tmp"
    if tmpfile.exists():
        print(f"Found temporary file {tmpfile}")
        second_tmpfile = outputdir/"thin_todel.tmp"
        if not second_tmpfile.exists():
            print("However, did not find second temp file 'thin_todel.txt'. Exiting!")
            exit(1)
        print("Continuing with last settings...")
        
        
        new_inds = np.genfromtxt(outputdir/"thin_newinds.tmp", dtype="int")
        inds_to_del = np.genfromtxt(outputdir/"thin_todel.tmp", dtype="int")
    else:
        new_inds = inds[::N]
        new_inds = np.append(new_inds, inds[-1])
        inds_to_del = [n for n in inds if not n in new_inds]

    print(f"Thinning down to {len(new_inds)} snapshots.")
    print("The following snapshots will be retained:")
    print(new_inds)
    print_size(outputdir, inds[0], len(new_inds))
    
    if not force:
        answer = input("Do you want to proceed, delete data, and adjust the list files?\n[type 'y' if yes]: ")
        if answer != "y":
            print("Abort.")
            exit(0)

    # print("Inds to be deleted:", inds_to_del)
    np.savetxt(outputdir/"thin_newinds.tmp", new_inds, fmt="%d")
    np.savetxt(outputdir/"thin_todel.tmp", inds_to_del, fmt="%d")
    
    print("Deleting snapshots...")
    for k, n in enumerate(inds_to_del):
        Ntodel = len(inds_to_del)
        print(f"\r{n} ({k+1} / {len(inds_to_del)}, {k/Ntodel*100:.1f}%)", end="", flush=True)
        p = outputdir / "snapshots" / f"{n}"
        if not p.exists():
            continue
        shutil.rmtree(p)
        
    print("\nWriting new snapshot list...")
    np.savetxt(outputdir/"snapshots"/"list.txt", new_inds, fmt="%d")
    
    print("Writing new time file...")
    lines = []
    with open(outputdir/"snapshots"/"timeSnapshot.dat", "r") as infile:
        for line in infile:
            if line.strip()[0] == "#":
                lines.append(line)
                continue
            ind = int(line.strip().split()[0])
            if ind in new_inds:
                lines.append(line)
    with open(outputdir/"snapshots"/"timeSnapshot.dat", "w") as outfile:
        outfile.writelines(lines)
    
    os.remove(outputdir/"thin_newinds.tmp")
    os.remove(outputdir/"thin_todel.tmp")
    print("Done")

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