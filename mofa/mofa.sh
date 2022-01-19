#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate diff

# Multi-omics data integration
cat ${BASEDIR}/../rna/sample.count.info ${BASEDIR}/../atac/sample.info  | cut -f 1-6 | sort | uniq -d > sample.info
if [ -f ${BASEDIR}/../atac/atac.combined.tss.counts.gz ]
then
    if [ -f ${BASEDIR}/../rna/gene.count ]
    then
	Rscript mofa.R sample.info ${BASEDIR}/../rna/gene.count ${BASEDIR}/../atac/atac.combined.tss.counts.gz
    fi
fi
