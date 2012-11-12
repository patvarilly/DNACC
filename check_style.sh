#!/bin/bash
# Run this script to make sure that all source file adhere to PEP 8 standards
# I spent a *whole* afternoon making sure that "pep8" produced no errors
# for the code.  If you make a change or a new example, you *MUST* make
# sure that it gets a clean bill from pep8

# Handle the easy cases automatically
# (adapted from https://gist.github.com/1903033)

# Removes whitespace chars from blank lines
pep8 -r --select=W293 -q --filename=*.py . | xargs sed -i 's/^[ \r\t]*$//'

# Removes trailing blank lines from files
pep8 -r --select=W391 -q --filename=*.py . | xargs sed -i -e :a -e '/^\n*$/{$d;N;ba' -e '}'

# Squashes consecutive blanks lines into one
pep8 -r --select=E303 -q --filename=*.py . | xargs sed -i '/./,/^$/!d'

find dnacc examples setup.py simple_dnacc -name '*.py' | xargs pep8 --ignore=E125 -r