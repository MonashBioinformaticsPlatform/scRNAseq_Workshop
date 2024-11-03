#!/bin/bash

# Create docs/workshop.R
Rscript scripts/make_R_script.R


Rscript -e 'bookdown::render_book("index.Rmd")'

