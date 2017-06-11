## Configuration file for rnaseq_pipeline.sh
##
## Place this script in a working directory and edit it accordingly.
##
## The default configuration assumes that the user unpacked the 
## chrX_data.tar.gz file in the current directory, so all the input
## files can be found in a ./chrX_data sub-directory

#if these programs are not in any PATH directories, please edit accordingly:
HISAT2=$(which hisat2)
STRINGTIE=$(which stringtie)
SAMTOOLS=$(which samtools)

TEMPLOC="./tmp" #this will be relative to the output directory