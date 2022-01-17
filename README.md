# T-ALL analysis scripts

This repository contains downstream analysis scripts for generating count matrices, differential expression analysis and differential peak calling. The fundamental [ATAC-Seq](https://github.com/tobiasrausch/ATACseq) and [RNA-Seq](https://github.com/tobiasrausch/RNAseq) processing pipelines need to be run beforehand because the following files are required:

- *.count: RNA-Seq raw gene counts
- *.fpkm: RNA-Seq FPKM for each gene
- *.peaks: Single-sample peak calls in BED format

