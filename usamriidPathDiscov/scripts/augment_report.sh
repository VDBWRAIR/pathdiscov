#!/bin/bash

sample="mysample"

outputdir=$1
sample=$2

d=$( dirname $( readlink -m $0 ) )

function relpath() {
    RET=$(python -c "from os.path import *; print relpath(realpath(\"$1\"));")
}

if cd $outputdir; then

	mkdir -p output/tmp

	cat quality_filter/logs/*.e > output/qstats.txt

	if [ -e iterative_blast_phylo_1/reports/contig.$sample.top.report.txt ]; then

		cat iterative_blast_phylo_1/reports/contig.$sample.top.report.txt | sed '1d' > output/tmp/tmp.txt
		cat ray2_assembly_1/ray2_assembly_1.fasta | awk '{if ($1 ~ />/) {printf $1"\t"} else {print length}}' | sed 's|>||' > ray2_assembly_1/contig_len.txt

		$d/tableX.pl --one 3 --two 1 output/tmp/tmp.txt ray2_assembly_1/contig_numreads.txt | awk -F"\t" '{ if (NF==23) {print $0"\t"0} else {print}}' > output/tmp/tmp2.txt
		$d/tableX.pl --one 3 --two 1 output/tmp/tmp2.txt ray2_assembly_1/contig_len.txt > output/tmp/tmp3.txt

		cat iterative_blast_phylo_1/reports/contig.$sample.top.report.txt | head -1 | tr "\n" "\t" > output/$sample.aug.report.txt
		echo -e "#reads_map2contig\tcontig_len(bp)" >> output/$sample.aug.report.txt
		cat output/tmp/tmp3.txt >> output/$sample.aug.report.txt

		if cd output; then
            echo "Making symlinks in $outputdir/output"
			ln -sf ../iterative_blast_phylo_1/reports/contig.$sample.top.report.txt
			ln -sf ../iterative_blast_phylo_1/reports/contig.$sample.top.phylo.txt 
			ln -sf ../iterative_blast_phylo_2/reports/R1.$sample.top.phylo.txt 
			ln -sf ../iterative_blast_phylo_2/reports/R2.$sample.top.phylo.txt 
            relpath ../iterative_blast_phylo_1/iterative_blast_phylo_1.contig
            ln -sf $RET
            relpath ../ray2_assembly_1/ray2_assembly_1.fasta
            ln -sf $RET
			#ln -sf $( readlink -m ../iterative_blast_phylo_1/iterative_blast_phylo_1.contig )
			#ln -sf $( readlink -m ../ray2_assembly_1/ray2_assembly_1.fasta )
		fi
    else
        echo "Missing iterative_blast_phylo_1/reports/contig.$sample.top.report.txt"
	fi

	rm -rf output/tmp
else
    echo "$outputdir missing"
fi
