#!/bin/bash

# replace ids w numbers and fastajoinlines

input=$2		# in
output=$3		# out

cat ${input} | fastajoinlines | awk 'BEGIN{id=1}{if ($0 ~ /^>/) {print ">c"id; print "c"id"\t"substr($0,2) > "contig.id"; id++;} else {print}}' > ${output}
