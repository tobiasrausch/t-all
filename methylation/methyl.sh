#!/bin/bash

# Convert to beta values
Rscript beta.R

# Subset to cimp probes
cat methylation.beta.tsv | grep -w -Ff <(cut -f 1 paper.table ) > methylation.beta.cimp.tsv

# Zip
gzip  methylation.beta.cimp.tsv 
gzip  methylation.beta.tsv
