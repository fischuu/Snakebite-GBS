#!/bin/bash

round() {
    printf "%.2f" "$1"
}

file1=$1
file2=$2
sample=$3
# Make sure that the provided numbers are actually numbers.
   if ! [[ $sample =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then >&2 echo "$sample is not a number"; exit 1; fi
  
extension1="${file1##*.}"
extension2="${file2##*.}"
filename1="${file1%.*}"
filename2="${file2%.*}"

fn1=$filename1"_"$sample".fastq"
fn2=$filename2"_"$sample".fastq"

if [ $extension1 == "gz" ]; then
  gunzip $file1;
  file1=$filename1;
  filename1="${file1%.*}"
  fn1=$filename1"_"$sample".fastq"
fi
if [ $extension2 == "gz" ]; then
  gunzip $file2;
  file2=$filename2;
  filename2="${file2%.*}"
  fn2=$filename2"_"$sample".fastq"
fi

lines=$(wc -l < $file1)
echo $lines
echo $sample

if (( $(awk 'BEGIN {print ("'$sample'" <= 1)}') )); then
  sample=$(awk 'BEGIN {printf("%.0f", "'$sample'" * "'$lines'")}')
fi

echo $sample

paste $file1 $file1 | \
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' | \
awk -v k=$sample 'BEGIN{srand(systime() + PROCINFO["pid"]); }{ s=x++<k?x- 1:int(rand()*x);
                  if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | \
awk -F"\t" -v file1=$fn1 -v file2=$fn2 '{print $1"\n"$3"\n"$5"\n"$7 > file1;\
                                         print $2"\n"$4"\n"$6"\n"$8 > file2}'
                                         
if [ $extension1 == "gz" ]; then
  gzip $fn1;
  gzip $file1;
fi
if [ $extension2 == "gz" ]; then
  gzip $fn2;
  gzip $file2;
fi
