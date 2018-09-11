#!/bin/bash
# combine the htseq counts into 1 file per pair. Need to look up pair information in the original file from ENA.
# This is run with no arguments but the folder name needs to be changed 
FOLDER='/bi/scratch/Summer/batch2'
FILE='/bi/scratch/Summer/ena_single_ended_all_info.txt'
#FILE='head_info.txt'
# the SRR column seems to be constant in the file
SRR_COLUMN=$(head -n 1 $FILE | awk '{for(i=1;i<=NF;i++){if($i~/^SRR/){print i}}}')

# column 1 contains the group name, use this to find the 2 samples
cut -f1 $FILE| uniq | while read -r line; 

	# find the pairs associated with the group name, -v to pass variable to awk
	do SRRs=$(awk -v group=$line -v col=$SRR_COLUMN '$1==group {print $(col)}' $FILE); 

	#get the 1st sample
	srr1=$(echo $SRRs| cut -f1 -d" "); 
	#echo $srr1; 
	#get the 2nd sample	
	srr2=$(echo $SRRs| cut -f2 -d" "); 
	#echo $srr2; 

	# check the ht-seq files exist
	if [ -f ${FOLDER}/$srr1*htseq* ];then  
	
		#locate the htseq-count file for the samples
		file1=$(find ${FOLDER}/$srr1*htseq*);
		#echo "file found";
		#echo $file1;
		
		if [ -f ${FOLDER}/$srr2*htseq* ];then  
	
			file2=$(find ${FOLDER}/$srr2*htseq*);
			#echo "file found";
			#echo $file2;
			
			# check the number of genes are the same
			if(($(wc -l $file1 |cut -f1 -d" ") == $(wc -l $file2 |cut -f1 -d" "))); 
				then 
					#echo "the number of rows match"; 
				
					R_input=${srr1}_${srr2}_counts.txt; 
					
					# create a header line
					echo -e "gene_id\t${srr1}\tgene_id\t${srr2}" > $R_input;

					#join the pair of files together, remove any lines that don't contain valid gene identifiers
					paste $file1 $file2 | awk '$1 ~ /ENSMUSG/ {print $0}' >> $R_input;

			# close if statement
			fi;		
		fi;		
	fi;
	
# close do statement
done; 

# one big piped command

#cut -f1 $FILE| uniq | while read -r line ; do OUT=$(awk -v group=$line -v col=$SRR_COLUMN '$1==group {print $(col)}' $FILE); srr1=$(echo $OUT| cut -f1 -d" "); echo $srr1; srr2=$(echo $OUT| cut -f2 -d" "); echo $srr2; file1=$(find /bi/scratch/Summer/data_mapping/$srr1*htseq*);file2=$(find /bi/scratch/Summer/data_mapping/$srr2*htseq*);echo $file1;echo $file2; if (($(wc -l $file1 |cut -f1 -d" ") == $(wc -l $file2 |cut -f1 -d" "))); then echo "the number of rows match"; R_input=${srr1}_${srr2}_counts.txt; echo -e "gene_id\t${srr1}\tgene_id\t${srr2}" > $R_input;