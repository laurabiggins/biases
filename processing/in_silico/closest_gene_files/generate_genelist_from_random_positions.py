# script to generate random positions and find the closest gene
# the only input file required is a gtf file - this should contain all the genes and their locations in the genome 


import random
from collections import defaultdict

import getopt
import sys


output_filename = "closest_genes.txt"
# number of genes to generate per chromosome
number_of_genes_per_chr = 20
gtf_file = "/bi/group/bioinf/Laura_B/bias_analysis/gtf_files/mouse/Mus_musculus.GRCm38.80.gtf"

options, remainder = getopt.getopt(sys.argv[1:], 'o:n:g:', ['output=', 
                                                         'number_of_genes=',
                                                         'gtf_file=',
                                                         ])
for opt, arg in options:
    if opt in ('-o', '--output'):
        output_filename = arg
    elif opt in ('-n', '--number_of_genes'):
        number_of_genes_per_chr = arg
    elif opt in ('-g', '--gtf_file'):
        gtf_file = arg

print 'OUTPUT          :', output_filename
print 'NUMBER OF GENES :', number_of_genes_per_chr
print 'GTF FILE        :', gtf_file
print 'REMAINING       :', remainder


# dictionary of random positions, keys are chr, values are lists of positions
random_positions = {}

with open("chr_list_noMT.txt") as f:
    genome_sizes = f.readlines()
    
# generating the random positions
for line in genome_sizes:

    row = line.split("\t")
    
    chr = row[0]
    chr_length = int(row[1])

    string = "length of chr %s is %s"
    print(string %(chr, chr_length))
    
    # create a list of random positions for each chr
    random_pos = []
    
    for x in range(number_of_genes_per_chr):
        random_pos.append(random.randint(1,chr_length))
    
    random_positions[chr] = random_pos
    print ("%s random positions generated for chr %s" % (len(random_pos), chr))
    

# load in gtf file
with open(gtf_file) as f:
    gtf_file = f.readlines()
    

# a dictionary of lists of tuples to store gene info for each chromosome
# keys are chr, values are gene info
gene_info = defaultdict(list)
    
    
for line in gtf_file:
    
    line = line.rstrip()
    if not line.startswith("#"):
    
        split_line=line.split("\t")
        
        # check whether the line contains a gene
        if (split_line[2] == "gene"):
            
            # parse the detail column that contains the name
            details = split_line[8].split(";")
            gene_id = details[0].replace("gene_id ", "")
            gene_id = gene_id.replace('"', "")
            
            gene_name = details[2].replace("gene_name ", "")
            gene_name = gene_name.replace('"', "")
            gene_name = gene_name.replace(' ', "")
            
            #print gene_name
           
            # load the gene info into a tuple and add to the list
            # in the gtf column 0 is chr, 3 is start, 4 is end
            chr = split_line[0]
            info = (gene_id, gene_name,chr,int(split_line[3]),int(split_line[4]),)
            
            gene_info[chr].append(info) 
            
            
# create a list of tuples to contain the closest genes
closest_genes = []


# go through each chr
for chr in random_positions:
    print ("chr is %s" %chr)
        
    for random_pos in random_positions[chr]:
    
        closest_pos = 1000000000
        closest_gene = ""
        
        for gene in gene_info[chr]:

            # first check whether the random position is within a gene
            if (gene[3] < random_pos and gene[4] > random_pos):
                            
                #print_str = "%s is within %s, start is %s, end is %s"
                #print (print_str % (random_pos, gene[4], gene[1], gene[2]))
                closest_pos=0
                closest_gene = gene
                
                # exit the for loop as we've found an overlapping gene
                break
                
            diff1 = abs(int(gene[3] - random_pos))
            diff2 = abs(int(gene[4] - random_pos))
            
            min_diff = min([diff1,diff2])
            
            # if smaller than the previous minimum, replace it  
            if(min_diff <= closest_pos):
                closest_pos = min_diff
                closest_gene = gene

        closest_genes.append(closest_gene)
        
header_line = "\t".join(["gene_id","gene_name","chromosome","start","end"])
# write out the closest genes
with open(output_filename, "w") as f:
    
    f.write(header_line)
    f.write("\n")
    f.write("\n".join(["\t".join([str(g) for g in closest_gene]) for closest_gene in closest_genes]))
