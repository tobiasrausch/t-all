SHELL := /bin/bash

# Targets
TARGETS = .conda .mamba .diff
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.mamba: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y -n base -c conda-forge mamba && touch .mamba

.diff: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n diff bioconductor-deseq2=1.34.0 bioconductor-vsn && touch .diff

clean:
	rm -f *~

distclean: clean
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/
