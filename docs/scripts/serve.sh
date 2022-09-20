#!/usr/bin/env bash

cd $(dirname $(realpath $0))/..
source .venv/bin/activate

if [[ ! -d .venv ]]; then
    bash scripts/setup_sphinx_venv.sh
fi

if [[ ! -d build/html ]]; then
    make html
fi

python scripts/run_livereload.py