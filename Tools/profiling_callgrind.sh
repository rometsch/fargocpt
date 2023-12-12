#!/usr/bin/env bash

echo "This tool runs a the callgrind tool of valgrind."
echo "The simulation will run much slower than usual (expect factor 20+ slower)."
echo "You need to call the binary directly without the wrapper."

FILEDIR="$(realpath $(dirname $0))"
FARGOCMD="$FILEDIR/../bin/fargocpt_exe"

# Check if at least two arguments are provided
if [ $# -lt 3 ]; then
  echo "usage: $(basename $0) Nthread <fargo command>"
  exit 1
fi

# Assign the first two arguments to variables
Nthread="$1"

# Shift the first two arguments
shift 1

# Print the first two arguments
echo "Using $Nproc MPI processes and $Nthread OMP threads"

export OMP_NUM_THREADS="$Nthread"
valgrind --tool=callgrind $Nproc "$FARGOCMD" "$@"



echo "Now view the results using e.g. 'kcachegrind $(ls | grep callgrind.out | sort -V | tail -n 1)'"