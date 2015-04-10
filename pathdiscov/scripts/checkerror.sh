#!/bin/bash

# $1 proj output directory
# $2 bool (1 verbose, 0 not)

# to run e.g., 
# checkerror.sh project/output/dir

# need one argument ! 
if [ $#	== 2 ];	then 
	
	res=$( readlink -m $1 ) 

	for i in ${res}/*; do 
 
		if [ -e $i/logs ]; then 

			j=$( basename $i);
			echo "***"${j}"***";

			for k in $i/logs*/*; do

				echo -n "."
#				echo -e "\t\t\t\t\t\t\t"$k
				myname=$( basename $k );
				
				# .o files
#				if [[ "$k" =~ o$ ]]; then
				if [ $( echo $k | perl -ne '{ if ($_ =~ m/\.o/) {print "1"} else {print "0"}}' ) == "1" ]; then

					if [ $( cat $k | grep "\[START\]" | wc -l ) == 1 ]; then
						if [ $( cat $k | grep "\[END\]" | wc -l ) != 1 ]; then
							echo
							echo -e "\t\t\t\t\t\t\t"$k
							echo "[error] job didnt finish: "${myname};
						fi;
					fi;
					
					if [ $( cat $k | grep "\[start\]" | wc -l ) != $( cat $k | grep "\[end\]" | wc -l ) ]; then
						echo
						echo -e "\t\t\t\t\t\t\t"$k
						echo "[error] number of start and end tags dont match: "${myname};
					fi;

					if [ "$2" == "1" ]; then echo "@@@ output ${j} > @@@"; echo; cat $k; echo "@@@ < output ${j} @@@"; echo; fi
					
				fi

				# .e files
#				if [[ "$k" =~ e$ ]]; then 
				if [ $( echo $k | perl -ne '{ if ($_ =~ m/\.e/) {print "1"} else {print "0"}}' ) == "1" ]; then

					ls -hl $k | awk -v name=${myname} -v k=${k} '$5!=0{printf "\n"; print "\t\t\t\t\t\t\t"k; print "[error] error file non-zero: "name}'
					
					if [ "$2" == "1" ]; then echo "@@@ error ${j} > @@@"; echo; cat $k; echo "@@@ < error ${j} @@@"; echo; fi
			
				fi

			done;
			
			echo
  
		fi; 
	done

else

	echo "error - need to supply an argument"

fi
