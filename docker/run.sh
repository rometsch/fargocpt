#!/usr/bin/env bash

COMMAND=$1
INPUTFILE="$(realpath $2)"
OUTPUTDIR="$(realpath $3)"

sudo docker run \
    -e LOCAL_USER_ID=$UID \
    -v $OUTPUTDIR:/output \
    -v $INPUTFILE:/app/setup.yml \
    fargocpt \
    fargocpt -np 1 \
    $COMMAND \
    setup.yml