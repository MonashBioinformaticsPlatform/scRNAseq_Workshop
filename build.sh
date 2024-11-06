#!/bin/bash

set -eu

# Start with a clean docs directory
rm -rf docs
mkdir docs

# Create docs/workshop.R
/persistent/software/apps/R/4.4.0/bin/Rscript scripts/make_R_script.R


/persistent/software/apps/R/4.4.0/bin/Rscript -e 'bookdown::render_book("index.Rmd")'

