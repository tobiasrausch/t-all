#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate diff

# Aggregate peaks
rm -f bam.lst
rm -f peak.lst
echo -e "id\tpatient\treltype\tsex\tfusion\tinirel\tfile" > sample.info
while read SAMPLE PSEUDO RELTYPE SEX FUSION
do
    if [ ${PSEUDO} == "P27" ]; then continue; fi
    for TYPE in INI REL
    do
	NUM=1
	for BAM in ../../bam/atac/PDX${SAMPLE}*${TYPE}*/*.final.bam
	do
	    PEAK=`echo ${BAM} | sed 's/.final.bam$/.peaks/'`
	    if [ -f ${PEAK} ]
	    then
		IDNAME=${PSEUDO}REP${NUM}${TYPE}PDX
		echo ${IDNAME}
		echo -e "${IDNAME}\t${PSEUDO}\t${RELTYPE}\t${SEX}\t${FUSION}\t${TYPE}\t${PEAK}" >> sample.info
		ln -s ${BAM} ${IDNAME}.bam
		ln -s ${BAM}.bai ${IDNAME}.bam.bai
		ln -s ${PEAK} ${IDNAME}.peaks
		echo ${IDNAME}.bam >> bam.lst
		echo ${IDNAME}.peaks >> peaks.lst
		NUM=`expr ${NUM} + 1`
	    fi
	done
    done
done < sample.mapping.pdx

# Create atac count matrix
/opt/dev/ATACseq/src/count.sh hg19 peaks.lst bam.lst atac.combined
rm -f *.bam *.bai *.peaks peaks.lst bam.lst

# T-ALL type 1 and type 2
Rscript diffATAC.R atac.combined.counts.gz sample.info all

# T-ALL type 1
Rscript diffATAC.R atac.combined.counts.gz sample.info 1

# T-ALL type 2
Rscript diffATAC.R atac.combined.counts.gz sample.info 2
