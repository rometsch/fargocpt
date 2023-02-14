#!/usr/bin/env python3
"""Thin out the output of a fargocpt simulation. E.g. keep only every 100th snapshot.
"""

import os
import argparse
import numpy as np
from pathlib import Path
import shutil
from typing import List


def main():
    opts = parse_args()

    if not hasattr(opts, "Nptrn"):
        mode = "inspect"
    elif hasattr(opts, "destination"):
        mode = "extract"
    else:
        mode = "delete"

    if mode == "inspect":
        inspect(opts.outputdir)
    elif mode == "extract":
        thin(opts.source, opts.Nptrn, force=opts.y,
             new_outputdir=opts.destination)
    elif mode == "delete":
        thin(opts.outputdir, opts.Nptrn, force=opts.y)


def thin(outputdir: Path, Nptrn: slice, force: bool = False, new_outputdir: Path = None):
    """Thin out an output directory or copy a reduced version.

    Delte all but every Nth snapshot or copy only scalar data and every Nth snapshot if new_outputdir is defined.

    Args:
        outputdir (Path): Path to the output dir.
        Nptrn (slice): Keep every Nth snapshot.
        force (bool, optional): Skip asking for confirmation. Defaults to False.
        new_outputdir (Path, optional): Path to the new directory where to copy data. Defaults to None.

    Raises:
        FileExistsError: if new_outputdir is not None, exists and is not empty.
    """
    outputdir = Path(outputdir)
    if new_outputdir is None:
        new_outputdir = outputdir
    else:
        new_outputdir = Path(new_outputdir)

        if new_outputdir.exists():
            if len(os.listdir(new_outputdir)) > 0:
                raise FileExistsError(
                    f"New outputdir {new_outputdir} already exists and is not empty!")
        else:
            os.makedirs(new_outputdir)

    inds = get_inds(outputdir)
    print(f"Output at path {outputdir}")

    tmpfile = new_outputdir/"thin_newinds.tmp"
    if tmpfile.exists():
        print(f"Found temporary file {tmpfile}")
        second_tmpfile = new_outputdir/"thin_todel.tmp"
        if not second_tmpfile.exists():
            print("However, did not find second temp file 'thin_todel.txt'. Exiting!")
            exit(1)
        print("Continuing with last settings...")

        new_inds = np.genfromtxt(outputdir/"thin_newinds.tmp", dtype="int")
        inds_to_del = np.genfromtxt(outputdir/"thin_todel.tmp", dtype="int")
    else:
        if not ":" in Nptrn:
            new_inds = [int(Nptrn)]
        else:
            try:
                istart, istop, istep = [
                    int(s) if s != "" else None for s in Nptrn.split(":")]
            except ValueError:
                istart, istop = [
                    int(s) if s != "" else None for s in Nptrn.split(":")]
                istep = 1
            sl = slice(istart, istop, istep)
            new_inds = inds[sl]
            if istop is None and inds[-1] not in new_inds:
                new_inds = np.append(new_inds, inds[-1])
        inds_to_del = [n for n in inds if not n in new_inds]

    print(f"Thinning down to {len(new_inds)} snapshots.")
    print("The following snapshots will be retained:")
    print(new_inds)
    print_size(outputdir, inds[0], len(new_inds))

    if not force:
        if outputdir == new_outputdir:
            query_str = "Do you want to proceed, delete data, and adjust the list files?\n[type 'y' if yes]: "
        else:
            query_str = "Do you want to proceed to copy the data?\n[type 'y' if yes]: "
        answer = input(query_str)
        if answer != "y":
            print("Abort.")
            exit(0)

    np.savetxt(new_outputdir/"thin_newinds.tmp", new_inds, fmt="%d")
    np.savetxt(new_outputdir/"thin_todel.tmp", inds_to_del, fmt="%d")

    if outputdir == new_outputdir:
        print("Deleting snapshots...")
        delete_snapshots(outputdir, inds_to_del)
    else:
        copy_output(outputdir, new_outputdir, new_inds)

    modify_time_files(outputdir, new_inds, new_outputdir=new_outputdir)

    os.remove(new_outputdir/"thin_newinds.tmp")
    os.remove(new_outputdir/"thin_todel.tmp")
    print("Done")


def copy_output(outputdir: Path, new_outputdir: Path, inds: List[int]):
    """Copy the content of an output dir to a new directory only keeping some snapshots.

    Args:
        outputdir (Path): Source output dir.
        new_outputdir (Path): Destination directory.
        inds (List[int]): Indices of snapshots to copy.
    """
    # copy everything but snapshots
    for p in os.listdir(outputdir):
        if p == "snapshots":
            continue
        src = outputdir / p
        dst = new_outputdir / p
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy(src, dst)
    # copy selected snapshots
    os.makedirs(new_outputdir / "snapshots")
    for n in inds:
        src = outputdir / "snapshots" / f"{n}"
        dst = new_outputdir / "snapshots" / f"{n}"
        shutil.copytree(src, dst)


def delete_snapshots(outputdir: Path, inds_to_del: List[int], print_progress: bool = True):
    """ Modify the snapshot list and time file to only include the retained snapshots.

    Args:
        outputdir (str): Path of the outputdir containing the directory 'snapshots'.
        inds_to_del (List[int]): List of indices to be deleted.
        print_progress (bool) [True]: Print progress of deletion.
    """
    for k, n in enumerate(inds_to_del):
        Ntodel = len(inds_to_del)
        if print_progress:
            print(
                f"\r{n} ({k+1} / {len(inds_to_del)}, {k/Ntodel*100:.1f}%)", end="", flush=True)
        p = outputdir / "snapshots" / f"{n}"
        if not p.exists():
            continue
        shutil.rmtree(p)


def modify_time_files(outputdir: Path, new_inds: List[int], new_outputdir: Path = None):
    """ Modify the snapshot list and time file to only include the retained snapshots.

    Args:
        outputdir (Path): Path of the outputdir containing the directory 'snapshots'.
        new_inds (List[int]): List of indices to keep.
        new_outputdir (Path, optional): Path to new outputdir. Defaults to None.
    """

    if new_outputdir is None:
        new_outputdir = outputdir

    print("\nWriting new snapshot list...")
    np.savetxt(new_outputdir/"snapshots"/"list.txt", new_inds, fmt="%d")

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
    with open(new_outputdir/"snapshots"/"timeSnapshot.dat", "w") as outfile:
        outfile.writelines(lines)


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

    thin_parser = subparsers.add_parser("delete")
    thin_parser.add_argument("outputdir", type=str)
    thin_parser.add_argument(
        "Nptrn", type=str, help="Pattern for inds. E.g. 1:10:3 for inds from 1 to 10 in steps of 3, or ::5 for every 5th. Same pattern as in fancy indexing.")
    thin_parser.add_argument("-y", action="store_true",
                             help="Answer yes to all questions.")

    extract_parser = subparsers.add_parser("extract")
    extract_parser.add_argument("source", type=str)
    extract_parser.add_argument("destination", type=str)
    extract_parser.add_argument(
        "Nptrn", type=str, help="Pattern for inds. E.g. 1:10:3 for inds from 1 to 10 in steps of 3, or ::5 for every 5th. Same pattern as in fancy indexing.")
    extract_parser.add_argument(
        "-y", action="store_true", help="Answer yes to all questions.")

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


if __name__ == "__main__":
    main()
