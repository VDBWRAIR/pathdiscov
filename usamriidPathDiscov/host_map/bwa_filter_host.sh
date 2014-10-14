#!/bin/bash

# bwa align to a reference; then extract from the original fastq file(s) the reads that didn't map
# can parallelize here

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
			
			# create sai's			
			cmd="bwa aln ${db} ${r1} > R1.paired.sai"
			echo "[cmd] "$cmd
			eval ${cmd}
			
			cmd="bwa aln ${db} ${r2} > R2.paired.sai"
			echo "[cmd] "$cmd
			eval ${cmd}
			
			cmd="bwa sampe ${db} R1.paired.sai R2.paired.sai ${r1} ${r2} ${opts} > paired.sam"
			echo "[cmd] "$cmd
			eval ${cmd}
			
			time2=$( date "+%s" )
			echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))				
			
			samtools view -bS paired.sam > paired.bam
			
			echo "[stats] paired stats"
			samtools flagstat paired.bam
			
			# get unmapped read IDs
			cat paired.sam | awk '$3=="*"' | cut -f1 | sort -u | awk '{print "@"$0}' > paired.unmap.id
			rm paired.sam
			
			# extract the 4 line chunks from a fastq file, given in argument 1 file, for IDs give in argument 2 file
			cmd="${d}/fastq_extract_id.pl ${r1} paired.unmap.id > R1.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}	

			cmd="${d}/fastq_extract_id.pl ${r2} paired.unmap.id > R2.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}		 		
		
		else # paired but not well-paired
		
			# find pairs and singletons		
			cmd="${d}/perlscripts_wrapper.pl get_common_uneven_files ${r1} ${r2} R1.single.fastq R1.paired.fastq R2.single.fastq R2.paired.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}	
			
			# paired -------------------------------			
			
			time1=$( date "+%s" ) 
			
			# create sai's			
			cmd="bwa aln ${db} R1.paired.fastq > R1.paired.sai"
			echo "[cmd] "$cmd
			eval ${cmd}
			
			cmd="bwa aln ${db} R2.paired.fastq > R2.paired.sai"
			echo "[cmd] "$cmd
			eval ${cmd}
			
			cmd="bwa sampe ${db} R1.paired.sai R2.paired.sai R1.paired.fastq R2.paired.fastq ${opts} > paired.sam"
			echo "[cmd] "$cmd
			eval ${cmd}
			
			time2=$( date "+%s" )
			echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))				
			
			samtools view -bS paired.sam > paired.bam
			
			echo "[stats] paired stats"
			samtools flagstat paired.bam
			
			# get unmapped read IDs
			cat paired.sam | awk '$3=="*"' | cut -f1 | sort -u | awk '{print "@"$0}' > paired.unmap.id
			rm paired.sam
			
			# extract the 4 line chunks from a fastq file, given in argument 1 file, for IDs give in argument 2 file
			cmd="${d}/fastq_extract_id.pl R1.paired.fastq paired.unmap.id > R1.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}	

			cmd="${d}/fastq_extract_id.pl R2.paired.fastq paired.unmap.id > R2.unmap.fastq"
			echo "[cmd] "$cmd
			eval ${cmd}	


			if [ -s R1.single.fastq ]; then
			
				time1=$( date "+%s" ) 
				
				# create sai's			
				cmd="bwa aln ${db} R1.single.fastq > R1.single.sai"
				echo "[cmd] "$cmd
				eval ${cmd}				
				
				cmd="bwa samse ${db} R1.single.sai R1.single.fastq ${opts} > R1.single.sam"
				echo "[cmd] "$cmd
				eval ${cmd}	
				
				time2=$( date "+%s" )
				echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))	
			
				samtools view -bS R1.single.sam > R1.single.bam
				
				echo "[stats] R1.single stats"
				samtools flagstat R1.single.bam
				
				cat R1.single.sam | awk '$3=="*"' | cut -f1 | sort -u | awk '{print "@"$0}' > R1.single.unmap.id
				rm R1.single.sam
				
				cmd="${d}/fastq_extract_id.pl R1.single.fastq R1.single.unmap.id >> R1.unmap.fastq"
				echo "[cmd] "$cmd
				eval ${cmd}					
				
			fi
						
			if [ -s R2.single.fastq ]; then
			
				time1=$( date "+%s" ) 
				
				# create sai's			
				cmd="bwa aln ${db} R2.single.fastq > R2.single.sai"
				echo "[cmd] "$cmd
				eval ${cmd}				
				
				cmd="bwa samse ${db} R2.single.sai R2.single.fastq ${opts} > R2.single.sam"
				echo "[cmd] "$cmd
				eval ${cmd}
				
				time2=$( date "+%s" )
				echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))		
			
				samtools view -bS R2.single.sam > R2.single.bam
				
				echo "[stats] R2.single stats"
				samtools flagstat R2.single.bam
				
				cat R2.single.sam | awk '$3=="*"' | cut -f1 | sort -u | awk '{print "@"$0}' > R2.single.unmap.id
				rm R2.single.sam
				
				cmd="${d}/fastq_extract_id.pl R2.single.fastq R2.single.unmap.id >> R2.unmap.fastq"
				echo "[cmd] "$cmd
				eval ${cmd}	
							
			fi
		
		fi	# not well paired
		
	else # single reads
	
		time1=$( date "+%s" ) 			

		# create sai's			
		cmd="bwa aln ${db} ${r1} > R1.single.sai"
		echo "[cmd] "$cmd
		eval ${cmd}				
		
		cmd="bwa samse ${db} R1.single.sai ${r1} ${opts} > R1.single.sam"
		echo "[cmd] "$cmd
		eval ${cmd}
		
		time2=$( date "+%s" )
		echo "[echo] mapping complete. deltat: "$(( $time2 - $time1 ))	

		samtools view -bS R1.single.sam > R1.single.bam
		
		echo "[stats] R1.single stats"
		samtools flagstat R1.single.bam
		
		cat R1.single.sam | awk '$3=="*"' | cut -f1 | sort -u | awk '{print "@"$0}' > R1.single.unmap.id
		rm R1.single.sam
		
		cmd="${d}/fastq_extract_id.pl ${r1} R1.single.unmap.id > R1.unmap.fastq"
		echo "[cmd] "$cmd
		eval ${cmd}	
								
	fi # single
	
fi # cd
