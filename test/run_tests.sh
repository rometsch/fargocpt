#!/usr/bin/env bash
# Run all automatic tests

FILEDIR="$(dirname $(realpath $0))"
cd "$FILEDIR"
for DIR in $(find . -mindepth 1 -maxdepth 1 -type d)
do
    if [[ -e "$DIR/run_auto_test.sh" ]];
    then
        ./$DIR/run_auto_test.sh
    else
        echo "WARNING: $DIR contains no automatic test!"
    fi
done