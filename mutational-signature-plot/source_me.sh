#!/bin/bash

R_PATH="/software/team113/dermatlas/R/R-4.2.2/bin"
export PATH="${R_PATH:?empty-path-variable}:${PATH}"

which R

export R_LIBS="/software/team113/dermatlas/R/R-4.2.2/lib/R/library:/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_mutational_signatures/renv/library/R-4.2/x86_64-pc-linux-gnu"