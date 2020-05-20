#!/bin/bash
read=$2
filename=$(basename -- "$1")
extension="${filename##*.}"
filename="${filename%.*.*.*}"
zcat $1 | awk '{{print (NR%4 == 1) ? "@'$filename'_" ++i "/'$read'": $0}}' > $filename.$2.fastq
