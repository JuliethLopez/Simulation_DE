#!/bin/bash

# Module bash
. /usr/share/modules/init/bash
module load HGI/softpack/users/hn4/seurat5/7.0

# R code
Rscript DESEQ2.R

