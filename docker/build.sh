#!/usr/bin/env bash
# Run this script from within the repository dir

cd $(dirname $(dirname $(realpath $0)))
sudo docker build -f docker/Dockerfile -t fargocpt .