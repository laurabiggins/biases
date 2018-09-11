#!/bin/bash 
FILE=$1

short_head=$(head -n 1 $FILE)
echo -e "gene_name\tgene_id\t${short_head}" > ${FILE}_with_gene_name.txt
cut -f1 ${FILE} |  ~/scripts/get_gene_name.sh | paste - <(tail -n +2 $FILE) >> ${FILE}_with_gene_name.txt
