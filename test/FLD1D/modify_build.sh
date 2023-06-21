#!/usr/bin/env bash

touch ../../src/fld.cpp
make -C ../../src CLI_OPTIONS="-DCONSTANT_FLD_FLUXLIMITER"