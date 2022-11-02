#!/usr/bin/env bash
rm ../../fargo
make -C ../../src -j
../../fargo start dust_drift.yml