#!/usr/bin/env bash

# This test is dedicated to correctly recognizing and conditionally processing .cram/.sam/.bam input alignment info

input=$1
fname=$(basename "$input")
ext="${fname##*.}"
nameM="${fname%.*}"
echo $ext
echo $nameM
output=$2
sf=$3
rep=$4
reference=$5
GTFFILE=$6

echo "+++++++++++++++++++++++++++++++"
echo "+++++++++++++++++++++++++++++++"
echo "${output}"
OUTDIR="${output}_F:${sf}_R:${rep}"

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

if [[ ! -f rnaseq_ballgown.R ]]; then
   echo "ERROR: R script rnaseq_ballgown.R not found in current directory!"
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
BALLGOWNLOC=./ballgown

LOGFILE=./run.log

for d in "$TEMPLOC" "$ALIGNLOC" "$BALLGOWNLOC" ; do
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
fi

#$SAMTOOLS index ${ALIGNLOC}/${sample}.bam
#$SAMTOOLS flagstat ${ALIGNLOC}/${sample}.bam

#echo "..removing intermediate files"
# rm ${TEMPLOC}/*
#rm ${TEMPLOC}/${sample}.unsorted.bam

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Assemble transcripts (StringTie)"
$STRINGTIE -p $NUMCPUS -G ${GTFFILE} -o ${ALIGNLOC}/${nameM}.gtf \
 -l ${sample} ${ALIGNLOC}/${nameM}.bam

## merge transcript file
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Merge all transcripts (StringTie)"
ls -1 ${ALIGNLOC}/*.gtf > ${ALIGNLOC}/mergelist.txt

$STRINGTIE --merge -p $NUMCPUS -G  ${GTFFILE} \
    -o ${BALLGOWNLOC}/stringtie_merged.gtf ${ALIGNLOC}/mergelist.txt

## estimate transcript abundance
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> Estimate abundance for each sample (StringTie)"

isample=$(basename ${input})
sample="${sample%%.*}"
dsample="${sample%_*}"
echo "+++++++++++++++++++++++++++++++++++"
echo "DSAMPLE HERE: " $dsample
echo "SAMPLE HERE: " $sample
echo "+++++++++++++++++++++++++++++++++++"
if [ ! -d ${BALLGOWNLOC}/${dsample} ]; then
   mkdir -p ${BALLGOWNLOC}/${dsample}
fi
$STRINGTIE -e -B -p $NUMCPUS -G ${BALLGOWNLOC}/stringtie_merged.gtf \
-o ${BALLGOWNLOC}/${dsample}/${dsample}.gtf ${ALIGNLOC}/${nameM}.bam

awk -F '\t|;|"| ' '$3 == "transcript" && match($19,/reference_id*/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11"\t"$16"\t"$21"\t"$26"\t"$31"\t"$36"\t"$41"\t"$46}' ${ALIGNLOC}/${nameM}.gtf > ${ALIGNLOC%/*}/transcripts.gtf

echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> DONE."
} #pipeline end

pipeline 2>&1 | tee $LOGFILE
