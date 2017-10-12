# An r script that runs gProfileR on a set of gene lists
# Takes gene info files as input and extracts the gene names to use in GO analysis

# Usage: Rscript run_gprofiler.r analysis_name files.genfo

# produces a histogram and results table of any GO categories that came up frequently

library("gProfileR")
library("ggplot2")

args <- commandArgs(trailingOnly = TRUE)

if (length(args)<=1) {
  stop("Files to be processed must be supplied", call.=FALSE)
} else if (length(args)>1) {
  
  analysis.name <- args[1]
   
  files <- args[2:length(args)]  
}

print(paste(length(files), "files to import"))
print(paste("starting", analysis.name, "analysis"))

# we want this to import a list of genes so we need a column named gene_name
first_line <- scan(files[1], nlines=1, what="character")


importFiles <- function(first_line, files){
  
  if(!sum(first_line == "gene_name")==1) stop ("couldn't find a column named 'gene_name'")
  
  if(length(first_line) == 1){
    gene.lists <- sapply(files, read.delim)
  }	
  else{
    
    col_no <- which(first_line == "gene_name")
    
    cc <- c(rep("NULL",col_no-1),"character",rep("NULL",length(first_line)-col_no))
    
    gene.lists <- sapply(files, read.delim, colClasses =cc)
  }
  message(paste("imported", length(gene.lists), "files"))
  gene.lists
}

gene.lists <- importFiles(first_line, files)

# run gProfileR, I don't want all the info that it returns so am selecting relevant columns
gprofiler.results <- lapply(gene.lists, function(x){
  gprofiler(as.vector(x), organism = "mmusculus", exclude_iea = TRUE, src_filter="GO", min_set_size=10, 
               max_set_size=2000)[,c(3:6,9,10,12)]
})


# all the significant GO categories that were returned for each genelist
go.ids <- unlist(sapply(gprofiler.results, "[", "term.id"))
go.terms <- unlist(sapply(gprofiler.results, "[", "term.name"))
category.size <- unlist(sapply(gprofiler.results, "[", "term.size"))

# count up how many times each category was identified
tabled.ids <- table(go.ids)
tabled.ids <- tabled.ids[order(tabled.ids, decreasing = TRUE)]

unique.locations <- !duplicated(go.ids)

unique.ids <- go.ids[unique.locations]
unique.name <- go.terms[unique.locations]
unique.size <- category.size[unique.locations]

id.order <- match(names(tabled.ids), unique.ids)
go.terms <- as.vector(unique.name[id.order])
go.size <- as.vector(unique.size[id.order])

df <- data.frame(id=names(tabled.ids), term=go.terms, freq=as.vector(tabled.ids), category.size=unique.size)

table.name <-  paste(analysis.name, "res_table.txt", sep="_")
write.table(x = df, file = table.name, quote = FALSE, sep="\t", row.names = FALSE)

if(nrow(df)>=1){

	# histogram of category frequency
	plot.name <- paste(analysis.name, "histogram.pdf", sep="_")
	title.text <- paste("Frequency of GO categories, no of gene lists =", length(gprofiler.results))
	pdf(file=plot.name)
	ggplot(df, aes(freq))+geom_histogram(binwidth=1)+labs(title=title.text)
	hist(tabled.ids, xlab="no of genelists in which GO category was overrepresented", main=title.text)
	dev.off()
}
