#!/usr/bin/env bash

cd source/Examples
# check is last command was successful otherwise exit
if [ $? -ne 0 ]; then
    echo "Error: Could not change directory to source/Examples"
    exit 1
fi

examples_dir="../../../examples"
# convert all ipynb files to markdown
for F in $(find "$examples_dir" -name "*.ipynb"); do
    jupyter nbconvert --to markdown --output-dir="." "$F"
done

# # write the markdown index
# echo "# Examples" > index.md
# echo "" >> index.md
# echo "This directory contains a collection of examples that illustrate the use of FargoCPT using Jupyter Notebooks." >> index.md
