#!/usr/bin/env bash
# Run this script from within the repository dir

cd $(dirname $(dirname $(realpath $0)))
docker build --build-arg STEP=false -f docker/Dockerfile -t fargocpt-dev .
