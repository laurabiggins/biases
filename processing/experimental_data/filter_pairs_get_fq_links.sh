#!/bin/bash
# This script processes sample names and information already downloaded from ena, selecting samples that have 2 replicates
# output is a list of locations for fastq files
FILE=$1
sed -i 's/ /_/g' $FILE

# select lines with rep1 and rep2, remove lines that contain single cell or tech rep 
# Some of the sample names are in different columns so loop through each field, if a match for rep1/2 is found print sample name (printf doesn't insert newline character), then print line, then break loop so we only get one match per line
# only keep entries with fastq files, sort by sample name, remove any entries with duplicate sample names
grep "rep1$\|rep2$" $FILE | grep -iv "singlecell\|single_cell\|single-cell\|techrep\|tech_rep\|tech-rep" | awk '{for(i=1;i<=NF;i++){ if($i~/\S*rep1$|\S*rep2$/){printf $i;printf "\t";print $0;break}}}' | awk '$0~/ftp/ && $0~/gz/' | sort - | awk '!seen[$1]++' > ${FILE}_unique_names.txt

# Add a new column with the sample names but removing the rep1/rep2 string so we can check for pairs
# sort by the first column, then the 3rd, not by the rep1/2 as matching the rep1/2 was beating the +-
awk '{print $1}' ${FILE}_unique_names.txt | sed 's/_rep1\|-rep1\|_rep2\|-rep2//' | sed 's/^GSM\w*://' | sed 's/^_//' | paste -  ${FILE}_unique_names.txt | sort -k 1,1 -k 3,3 - >  ${FILE}_sorted.txt

# only keep the pairs
cut -f1 ${FILE}_sorted.txt | uniq -c | awk '$1==2 {print $2}' | awk 'FNR==NR{a[$1];next}($1 in a){print}' - ${FILE}_sorted.txt > ${FILE}_all_info.txt

# extract the ftp links - some of the entries have the file locations in different columns so we check each column and break the loop once found so we only get one file per sample
# prepend ftp:// to each line so that clusterflow can process the files
awk '{for(i=1;i<=NF;i++) {if($i ~ /^ftp/ && $i ~ /gz/){print $i; break}}}' ${FILE}_all_info.txt | awk '{print "ftp://" $0}' > ${FILE}_fq_links.txt

# clean up
rm ${FILE}_unique_names.txt ${FILE}_sorted.txt
rename ".txt_" "_" ${FILE}_all_info.txt ${FILE}_fq_links.txt
