#!/bin/bash

input=$1
flanking=$2
output1=$3
output2=$4
inputvcf=$5

# Prepare the loci file
sed '/^#/ d' $input | cut -f1,2 | awk -v s=$flanking -v OFS= '{print $1,":", $2-s,"-",$2+s}' > $output1
samtools faidx $inputvcf -n 300 -r $output1 > $output2
