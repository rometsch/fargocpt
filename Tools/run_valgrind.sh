#!/usr/bin/env bash
SETUP="$1"
valgrind \
    --leak-check=full \
    --show-leak-kinds=all \
    --track-origins=yes \
    ./fargo start $SETUP \
    2>&1 | tee valgrind.log