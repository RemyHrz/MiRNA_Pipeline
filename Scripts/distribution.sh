#!/bin/bash

filename=$(basename "$1" .sam)
out=$3
common_name=$2
filetmp="${out}.tmp"


samtools view -@ $3 -F 4 $1 | cut -f 1,10 | awk '{FS=OFS="\t";print $1,length($2)}' | sort | uniq > $filetmp

count=$(wc -l $filetmp|cut -d" " -f 1)

cut -f 2 $filetmp | sort | uniq -c |awk '{print $2, $1}' > $out



