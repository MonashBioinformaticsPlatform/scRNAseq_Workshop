#!/bin/bash

set -eu

# Start with a clean docs directory
rm -rf docs
mkdir docs

# Create docs/workshop.R
Rscript scripts/make_R_script.R


Rscript -e 'bookdown::render_book("index.Rmd")'

