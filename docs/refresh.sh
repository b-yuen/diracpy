#!/bin/bash

# Run the make html command twice
make clean
make html

# Open the resulting HTML in Safari
open -a safari build/html/index.html
