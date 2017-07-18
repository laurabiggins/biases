#!/bin/bash
# This script takes a mapped bam file, sorts the bam file, uses htseq-count to count reads over genes, then appends gene symbols to the output.
# to run this  
# qsub -cwd -V -l h_vmem=10G ./run_htseq.sh mapped_file.bam
for FILE in "$@"

do	
	
	echo "trying to sort and get counts for file $FILE"  > log.txt
	HT_OUTPUT='_htseq_output.txt'
	SORTED='sorted_'
	#echo "sorted file should be called $SORTED$FILE"  >> log.txt
	
	# samtools doesn't seem to get rid of the temp files unless it creates a merged version so I'm writing the file out. We can clean up afterwards.
	samtools sort -O BAM -o $SORTED$FILE $FILE
	
	# get the counts
	htseq-count --format bam --stranded=reverse $SORTED$FILE /bi/group/bioinf/Laura_B/bias_analysis/GRCm38_loci/Mus_musculus.GRCm38.88.gtf.gz > $FILE$HT_OUTPUT;
	
	# add in the gene names
	perl /bi/group/bioinf/Laura_B/bias_analysis/TIDIED/scripts/get_genesymbols_from_ensembl_ids.pl --gtf /bi/group/bioinf/Laura_B/bias_analysis/GRCm38_loci/Mus_musculus.GRCm38.88.gtf.gz --file $FILE$HT_OUTPUT;
	
done

