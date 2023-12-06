#!/usr/bin/env bash
make html
rsync -r build/html/ ../docs