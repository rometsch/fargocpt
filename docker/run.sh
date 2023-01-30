#!/usr/bin/env bash

if [ $# -ne 3 ]; then
    echo "usage: $0 {start/auto} {inputfile} {simulationdir}"
    exit 1
fi

COMMAND=$1
INPUTFILE="$(realpath $2)"
SIMDIR="$(realpath $3)"

if [ ! -e "$INPUTFILE" ]; then
    echo "Input file $INPUTFILE does not exist!"
    exit 1
fi

if [ ! -d "$SIMDIR" ]; then
    echo "Simulation directory $SIMDIR does not exist!"
    exit 1
fi

sudo docker run \
    -e LOCAL_USER_ID=$UID \
    -v $SIMDIR:/simulation \
    -v $INPUTFILE:/simulation/setup.yml \
    fargocpt \
    fargocpt -np 1 \
    $COMMAND \
    setup.yml