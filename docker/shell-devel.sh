#!/usr/bin/env bash

if (("$#"<"1")); then
    echo "usage: $0 {simulationdir}"
    exit 1
fi

SIMDIR="$(realpath $1)"

if [ ! -d "$SIMDIR" ]; then
    echo "Simulation directory $SIMDIR does not exist!"
    exit 1
fi

docker run \
    --user $(id -u):$(id -g) \
    -v $SIMDIR:/project \
    --rm \
    -it \
    fargocpt-devel \
    /bin/bash
