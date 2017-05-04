#!/usr/bin/env bash

# This test is dedicated to correctly recognizing and conditionally processing .cram/.sam/.bam input alignment info

input=$1
output=$2
sf=$3
rep=$4
reference=$5
GTFFILE=$6
NUMCPUS=$7

fname=$(basename "$input")
ext="${fname##*.}"
nameM="${fname%.*}"

if [[ ${sf#*.} -eq 0 ]]; then
    sfName="1.${sf#*.}"
else
    sfName="0.${sf#*.}"
fi

if [[ $sfName == "0.0" ]]; then
    sfName="1.0"
fi
OUTDIR="${output}_F:${sfName}_R:${rep}"

## load variables
if [[ ! -f ./rnaseq_pipeline.config.sh ]]; then
 usage
 echo "Error: configuration file (rnaseq_pipeline.config.sh) missing!"
 exit 1
fi

source ./rnaseq_pipeline.config.sh
WRKDIR=$(pwd -P)
errprog=""
if [[ ! -x $SAMTOOLS ]]; then
    errprog="samtools"
fi
if [[ ! -x $HISAT2 ]]; then
    errprog="hisat2"
fi
if [[ ! -x $STRINGTIE ]]; then
    errprog="stringtie"
fi

if [[ "$errprog" ]]; then    
  echo "ERROR: $errprog program not found, please edit the configuration script."
  exit 1
fi

#determine samtools version
newsamtools=$( ($SAMTOOLS 2>&1) | grep 'Version: 1\.')

set -e
#set -x

if [[ $OUTDIR != "." ]]; then
  mkdir -p $OUTDIR
  cd $OUTDIR
fi

SCRIPTARGS="$@"
ALIGNLOC=./hisat2

LOGFILE=./run.log

for d in "$TEMPLOC" "$ALIGNLOC" ; do
 if [ ! -d $d ]; then
    mkdir -p $d
 fi
done

# main script block
pipeline() {

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 $SCRIPTARGS

isample=$(basename ${input})
sample=$input

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Alignments conversion (SAMTools)"
if [[ "$newsamtools" ]]; then
 
 if [[ "$ext" = "sam" ]]; then
    echo "sam"
    $SAMTOOLS view -h -s ${sf} -S -b ${TEMPLOC}/${nameM}.sam | \
     $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${nameM}.bam -

 elif [[ "$ext" = "bam" ]]; then
    echo "bam"
 elif [[ "$ext" = "cram" ]]; then
    echo "cram"
    if [ ! -f ${output}.sam ]; then
        $SAMTOOLS view -h -T ${reference} -o ${output}.sam ${input}
    fi
    $SAMTOOLS view -h -s ${sf} -S -b ${output}.sam | \
     $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${nameM}.bam -
else
    echo "unknown fileformat. please provide .sam/.bam/.cram"
 fi
   
else
 if [[ "$ext" = "sam" ]]; then
    echo "sam"
    if [ ! -f ${output}.sam ]; then
        cp ${input} ${output}.sam
    fi

    $SAMTOOLS view -h -s ${sf} -S -b ${output}.sam | \
     $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${nameM}.bam -

 elif [[ "$ext" = "bam" ]]; then
    echo "bam"
 elif [[ "$ext" = "cram" ]]; then
    echo "cram"
    if [ ! -f ${output}.sam ]; then
        $SAMTOOLS view -h -T ${reference} -o ${output}.sam ${input}
    fi
    $SAMTOOLS view -h -s ${sf} -S -b ${output}.sam | \
     $SAMTOOLS sort -@ $NUMCPUS -o ${ALIGNLOC}/${nameM}.bam -
else
    echo "unknown fileformat. please provide .sam/.bam/.cram"
 fi
fi

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
$STRINGTIE -p $NUMCPUS -G ${GTFFILE} -e -o ${ALIGNLOC}/${nameM}.gtf \
 -l ${sample} ${ALIGNLOC}/${nameM}.bam

awk -F '\t|;|"| ' '$3 == "transcript" && match($19,/ref_gene_name*/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11"\t"$16"\t"$21"\t"$26"\t"$31"\t"$36"\t"$41"\t"$46}' ${ALIGNLOC}/${nameM}.gtf > ${ALIGNLOC%/*}/transcripts.gtf

rm -rf ${ALIGNLOC}/*.bam

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
