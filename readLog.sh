#!/usr/bin/env bash

inputDir=$1
output=$2

touch $output
> $output

for dir in $inputDir*/ ; do
    dirName=$(basename "$dir")
    for file in $dir/*.log ; do
        sampleBase=$(basename "$file")
        sample="${sampleBase%.*}"
        nReads=$(head -n 1 $file | awk -F " " '{print $2}' -)
        percentAlignBase=$(tail -n 1 $file | awk -F " " '{print $1}')
        percentAlign="${percentAlignBase%\%*}"
        echo $dirName,$sample,$nReads,$percentAlign >> $output
    done
done
