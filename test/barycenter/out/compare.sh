#!/usr/bin/env bash
COMP="../../binary_io/print_differences.py"
if [[ "$#" == "0" ]]; then
	N="100"
else
	N="$1"
fi
TOCOMPARE="gasdens$N.dat"
$COMP master/$TOCOMPARE refractored/$TOCOMPARE
