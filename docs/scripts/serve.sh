#!/usr/bin/env bash

cd $(dirname $(realpath $0))/..

if [[ ! -d .venv ]]; then
    bash scripts/setup_sphinx_venv.sh
fi

if [[ ! -d build/html ]]; then
    make html
fi

source .venv/bin/activate
python scripts/run_livereload.py