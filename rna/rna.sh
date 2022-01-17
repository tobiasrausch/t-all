#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate diff

# Aggregate FPKM and raw counts
for FTYPE in fpkm count
do
    echo -e "id\tpatient\treltype\tsex\tfusion\tinirel\tfile" > sample.${FTYPE}.info
    zcat gene.table.gz > gene.${FTYPE}
    while read SAMPLE PSEUDO RELTYPE SEX FUSION
    do
	for TYPE in INI REL
	do
	    for GCOUNT in ../../bam/rna/${SAMPLE}*${TYPE}*/*.gene.${FTYPE}
	    do
		if [ -f ${GCOUNT} ]
		then
		    IDNAME=`echo ${GCOUNT} | sed 's/\/[^\/]*$//' | sed 's/^.*\///' | sed "s/^${SAMPLE}/${PSEUDO}/" | sed "s/^${PSEUDO}${TYPE}PDX/${PSEUDO}REP1${TYPE}PDX/"`
		    if [ `echo ${IDNAME} | grep -c "SECONDARY"` -eq 1 ]; then continue; fi
		    echo ${IDNAME}
		    echo -e "${IDNAME}\t${PSEUDO}\t${RELTYPE}\t${SEX}\t${FUSION}\t${TYPE}\t${GCOUNT}" >> sample.${FTYPE}.info
		    echo -e "gene\t${IDNAME}" > local.count
		    tail -n +2 ${GCOUNT} >> local.count
		    if [ -f gene.${FTYPE} ]
		    then
			diff <(cut -f 1 gene.${FTYPE}) <(cut -f 1 local.count)
			paste gene.${FTYPE} <(cut -f 2 local.count) > gene.${FTYPE}.tmp
			mv gene.${FTYPE}.tmp gene.${FTYPE}
		    fi
		    rm local.count
		fi
	    done
	done
    done < sample.mapping.pdx
done

# Differential expression REL vs. INI (all T-ALL types)
Rscript diffRNA.R gene.count sample.count.info all

# Differential expression REL vs. INI (T-ALL type 1)
Rscript diffRNA.R gene.count sample.count.info 1

# Differential expression REL vs. INI (T-ALL type 2)
Rscript diffRNA.R gene.count sample.count.info 2
