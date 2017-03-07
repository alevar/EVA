#!/usr/bin/env bash

#======================================
# Parameters previously in config.sh
#======================================
NUMCPUS=8
STRINGTIE=$(which stringtie)
SAMTOOLS=$(which samtools)

BASEDIR=$(pwd -P)/assembly

#FASTQLOC="$BASEDIR/samples"
#GENOMEIDX="$BASEDIR/indexes/chrX_tran"
#PHENODATA="$BASEDIR/geuvadis_phenodata.csv"

## list of samples 
## (only paired reads, must follow _1.*/_2.* file naming convention)
# reads1=(${FASTQLOC}/*_1.*)
# reads1=("${reads1[@]##*/}")
# reads2=("${reads1[@]/_1./_2.}")

errprog=""
if [[ ! -x $SAMTOOLS ]]; then
    errprog="samtools"
fi
if [[ ! -x $STRINGTIE ]]; then
    errprog="stringtie"
fi

if [[ "$errprog" ]]; then    
  echo "ERROR: $errprog program not found, please edit the parameters"
  exit 1
fi

#determine samtools version
newsamtools=$( ($SAMTOOLS 2>&1) | grep 'Version: 1\.')

set -e
#set -x

SCRIPTARGS="$@"
ALIGNLOC=./data

LOGFILE=./assembly/run.log

# main script block
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS

# Assemble the transcripts

sampleID=SRR821573 #sample ID

refann=/scratch0/genomes/hg38/annotation/hg38_p8.biotype_flt.cls.gff3

refseq=/scratch0/genomes/hg38/hg38.fa

samtools view -h -T $refseq ./data/$sampleID.cram | stringtie -p8 -m 150 -G $refann \
 -o ./assembly/$sampleID.gtf -A ./assembly/genes.tab -


# echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
# $STRINGTIE -p $NUMCPUS -G ${ALIGNLOC}/hg38_ucsc.annotated.gtf -o ./assembly/$sampleID.gtf \
#  -l SRR821573 ${ALIGNLOC}/SRR821573.cram

## merge transcript file
# echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Merge all transcripts (StringTie)"
# ls -1 ./assembly/*.gtf > ./assembly/mergelist.txt

# $STRINGTIE --merge -p $NUMCPUS -G  ${ALIGNLOC}/hg38_ucsc.annotated.gtf \
#     -o ./assembly/stringtie_merged.gtf ./assembly/mergelist.txt

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
