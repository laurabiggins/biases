#!/bin/bash
# This script takes a mapped bam file, sorts the bam file, uses htseq-count to count reads over genes
# to run this  
# for i in *bam; do echo $i; qsub -cwd -V -l h_vmem=16G -N htseq sh ../htseq_wrapper.sh $i;done
FILE=$1
GTF='/bi/scratch/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38.90.gtf'
samtools sort -O BAM $FILE | htseq-count --quiet --type gene --format bam --stranded=no - $GTF > ${FILE}_htseq_output.txt 
