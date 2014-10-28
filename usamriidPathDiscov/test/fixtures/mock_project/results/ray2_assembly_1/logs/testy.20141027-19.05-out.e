cat: /home/AMED/tyghe.vallard/Projects/usamriidPathDiscov/testy/results/ray2_assembly_1/results/Contigs.fasta: No such file or directory
Warning: Empty input file
Error: No unambiguous stretches of characters in the input.  Aborting...
Error: Encountered internal Bowtie 2 exception (#1)
Command: bowtie2-build /home/AMED/tyghe.vallard/Projects/usamriidPathDiscov/testy/results/ray2_assembly_1/ray2_assembly_1.fasta bowtie2_index/contigs 
Deleting "bowtie2_index/contigs.3.bt2" file written during aborted indexing attempt.
Deleting "bowtie2_index/contigs.4.bt2" file written during aborted indexing attempt.
Could not locate a Bowtie index corresponding to basename "bowtie2_index/contigs"
Error: Encountered internal Bowtie 2 exception (#1)
Command: /home/AMED/tyghe.vallard/Projects/usamriidPathDiscov/usamriidPathDiscov/bin/bowtie2-align --wrapper basic-0 -q -x bowtie2_index/contigs -S bowtie2_mapping/out.sam --local -1 R1.paired.fastq -2 R2.paired.fastq -U R1.single.fastq,R2.single.fastq 
bowtie2-align exited with value 1
