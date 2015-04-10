#!/bin/bash

# replace ids w numbers and fastajoinlines

d=$1			# scripts
input=$2		# in
output=$3		# out

cat ${input} | ${d}/fastajoinlines | awk 'BEGIN{id=1}{if ($0 ~ /^>/) {print ">c"id; print "c"id"\t"substr($0,2) > "contig.id"; id++;} else {print}}' > ${output}
