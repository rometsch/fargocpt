#!/usr/bin/env bash

if [ $# -le 2 ]; then
    echo "usage: $0 {simulationdir} {fargo args}"
    exit 1
fi

SIMDIR="$(realpath $1)"

if [ ! -d "$SIMDIR" ]; then
    echo "Simulation directory $SIMDIR does not exist!"
    exit 1
fi

shift 1

docker run \
    -e LOCAL_USER_ID=$UID \
    -v $SIMDIR:/project \
    fargocpt \
    fargocpt -np 1 "$@"
