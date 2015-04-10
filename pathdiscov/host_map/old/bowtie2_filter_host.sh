#!/bin/bash

# bowtie2 align to a reference; then extract from the original fastq file(s) the reads that didn't map

outputdir=$1			# outputdir
d=$2					# scripts dir
r1=$3					# reads 1 file
r2=$4					# reads 2 file
db=$5					# database
paired=$6				# paired - 1 if paired data, 0 if not 
wellpaired=$7			# wellpaired - 1 if paired data "well-paired" (i.e., mate files have exact same IDs), 0 if not
opts=$8					# options

if ! cd ${outputdir}; then
 
	echo "[error] changing directories"; 

else
	
	# paired data
	if [ "$paired" == 1 ]; then 

		# well paired 
		if [ "$wellpaired" == 1 ]; then 
		
			time1=$( date "+%s" )
			
			cmd="bowtie2 -q -x ${db} -1 ${r1} -2 ${r2} -S out.sam ${opts}"
			echo "[cmd] "$cmd
			eval ${cmd}
			
			time2=$( date "+%s" )
			echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))		
			
			samtools view -bS out.sam > out.bam
			
			echo "[stats] paired stats"
			samtools flagstat out.bam
			
			# get unmapped read IDs
			cat out.sam | awk '$3=="*"' | cut -f1 > tmp1.paired.unmap.id
			sort -u tmp1.paired.unmap.id > tmp2.paired.unmap.id
			cat tmp2.paired.unmap.id | awk '{print "@"$0}' > paired.unmap.id
			
			rm tmp1.paired.unmap.id tmp2.paired.unmap.id
			rm out.sam
			
			# extract the 4 line chunks from a fastq file, given in argument 1 file, for IDs give in argument 2 file		
			cmd="${d}/fastq_extract_id.pl ${r1} paired.unmap.id > R1.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}
					
			cmd="${d}/fastq_extract_id.pl ${r2} paired.unmap.id > R2.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}
	
			# these will automatically be well paired because of the way we got unmapped read IDs (i.e., if one maps and one doesnt it still goes into the unmap bin)
			
			# find pairs and singletons
			# input_R1 input_R2 R1.single.fastq R1.paired.fastq R2.single.fastq R2.paired.fastq
	#		cmd="${d}/perlscripts_wrapper.pl get_common_uneven_files R1.unmap.fastq R2.unmap.fastq R1.unmap.single.fastq R1.unmap.paired.fastq R2.unmap.single.fastq R2.unmap.paired.fastq"
	#		echo "[cmd] "$cmd
	#		eval ${cmd}		 		
		
		else # paired but not well-paired
		
			cmd="${d}/perlscripts_wrapper.pl get_common_uneven_files ${r1} ${r2} R1.single.fastq R1.paired.fastq R2.single.fastq R2.paired.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}		
			
			# paired & unpaired -------------------------------
			
			if [ -s R1.single.fastq -a -s R2.single.fastq  ]; then
				cmd="bowtie2 -q -x ${db} -1 R1.paired.fastq -2 R2.paired.fastq -U R1.single.fastq,R2.single.fastq -S out.sam ${opts}"
			elif [ -s R1.single.fastq ]; then
				cmd="bowtie2 -q -x ${db} -1 R1.paired.fastq -2 R2.paired.fastq -U R1.single.fastq -S out.sam ${opts}"
			elif [ -s R2.single.fastq ]; then
				cmd="bowtie2 -q -x ${db} -1 R1.paired.fastq -2 R2.paired.fastq -U R2.single.fastq -S out.sam ${opts}"
			else
				cmd="bowtie2 -q -x ${db} -1 R1.paired.fastq -2 R2.paired.fastq -S out.sam ${opts}"		
			fi
	
			time1=$( date "+%s" )
			
			echo "[cmd] "$cmd
			eval ${cmd}
			
			time2=$( date "+%s" )
			echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))
	
			samtools view -bS out.sam > out.bam
	
			echo "[stats]"
			samtools flagstat out.bam
						
			# get unmapped read IDs
			cat out.sam | awk '$3=="*"' | cut -f1 > tmp1.unmap.id
			sort -u tmp1.unmap.id > tmp2.unmap.id
			cat tmp2.unmap.id | awk '{print "@"$0}' > unmap.id
			
			rm tmp1.unmap.id tmp2.unmap.id
			rm out.sam			
			
			# extract the 4 line chunks from a fastq file, given in argument 1 file, for IDs give in argument 2 file
			cmd="${d}/fastq_extract_id.pl ${r1} unmap.id > R1.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}			
	
			cmd="${d}/fastq_extract_id.pl ${r2} unmap.id > R2.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}
	
		fi	# not well paired
	else # single reads
	
		time1=$( date "+%s" )

		cmd="bowtie2 -q -x ${db} -U ${r1} -S out.sam ${opts}"
		echo "[cmd] "$cmd
		eval ${cmd}
		
		time2=$( date "+%s" )
		echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))		
		
		samtools view -bS out.sam > out.bam
		
		echo "[stats]"
		samtools flagstat out.bam
		
		# get unmapped read IDs
		cat out.sam | awk '$3=="*"' | cut -f1 > tmp1.single.unmap.id
		sort -u tmp1.single.unmap.id > tmp2.single.unmap.id
		cat tmp2.single.unmap.id | awk '{print "@"$0}' > single.unmap.id
		
		rm tmp1.single.unmap.id tmp2.single.unmap.id
		rm out.sam
		
		cmd="${d}/fastq_extract_id.pl ${r1} single.unmap.id > R1.unmap.fastq"
		echo "[cmd] "$cmd
		eval ${cmd}							
	fi # single
fi # cd

