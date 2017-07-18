Scripts for processing public experimental data

------------------------------------
 mapping the fastq files - hisat2
------------------------------------

Using ENCODE rna-seq data (from a project already in Labrador - Encode_CSHL_2015). 

The data was mapped with hisat2 using clusterflow.

Usage:
cf --genome GRCm38 hisat2  ../FastQ/*gz

-------------------------------------------
 getting counts over genes - run_htseq.sh 
------------------------------------------

The bash script run_htseq.sh sorts the bam file using samtools sort, and passes the sorted bam file to htseq-count to count the reads over genes. 
htseq-count requires a gtf file to define the genes - download from ensembl http://www.ensembl.org/info/data/ftp/index.html
The version used here was Mus_musculus.GRCm38.88.gtf.gz.
A perl script - get_genesymbols_from_ensembl_ids.pl - then uses the Ensembl gene ID to look up the official gene symbol and adds this to the file.

Usage:
qsub -cwd -V -l h_vmem=10G run_htseq.sh mapped_file.bam

