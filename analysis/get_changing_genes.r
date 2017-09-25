# Rscript /bi/group/bioinf/Laura_B/bias_analysis/TIDIED/scripts/get_changing_genes.r  Encode_CSHL_sample_sheet.txt encode_rna_seq_analysis/*gene_names.txt

# sample sheet must have headers as below
# sample \t group \t description

# the data files must have the structure
# gene_name \t score \t ensembl_id

# this script calculates z-scores from local distributions and selects the top n (default 200) changing genes based on the z-scores


#rm(list=ls())

library("beanplot")

#setwd("D:/projects/biases/gene_count_files/")
#sample_sheet_file <- "../sample_sheet_2reps.txt"
#files <- list.files(pattern="SRR*")


args <- commandArgs(trailingOnly = TRUE)

#print(args[1])

#print ("=================")

if (length(args)<=1) {
  stop("Sample sheet and files to be processed must be supplied", call_=FALSE)
} else if (length(args)>1) {
  
  sample_sheet_file <- args[1]
  
  files <- args[2:length(args)]  
}

# read in the sample sheet
sample_sheet <- read.delim(sample_sheet_file)

# check which files match the sample sheet
matched_files <- unlist(sapply(sample_sheet$sample, function(x) grep(x, files, value=TRUE)))

matched_sample_names <- unlist(sapply(sample_sheet$sample, function(x){
  match <- grep(x, files, value=FALSE)
  ifelse(match>=1, return(as.vector(x)), return (NULL))
}))


#==========================
# sorting the sample sheet
#==========================

# remove any files not found from the sample sheet
sample_sheet <- sample_sheet[as.character(sample_sheet$sample) %in% matched_sample_names,]

# split by group
samples <- split(sample_sheet$sample, f = sample_sheet$group, drop=TRUE)

# remove the levels so they don't confuse things later on
samples <- lapply(samples, as.vector)


#=====================
# Importing the files
#=====================

# get our file names in the right structure
files_to_import <- lapply(samples, function(groups) sapply(groups, function(x) grep(x, files, value=TRUE)))

# import datasets
datasets <- lapply(files_to_import, function(x) lapply(x, read.delim))


# Collapse each group to a dataframe
# each group should have all the same gene names/ids so we should be able to collapse these to a dataframe
df <- lapply(datasets, function(x){
  
  sum_of_mismatches <- 0  
  
  # The ensembl ids should be the same for all the files
  for (i in 1:length(x))  for(j in 2:length(x)){
    sum_of_mismatches <- sum_of_mismatches + sum(x[[i]]$ensembl_id != x[[j]]$ensembl_id)    
  }
  if(sum_of_mismatches > 0){
    stop("ensembl ids don't match")
  }
  data.frame(ensembl_id = x[[1]]$ensembl_id, gene_name = x[[1]]$gene_name, sapply(x, `[[`, 'score'))
})


# convert to matrices for easier stats
mat_raw <- lapply(df, function(x){
  the_matrix <- as.matrix(x[,3:ncol(x)])
  row.names(the_matrix) <- (x[,"ensembl_id"])
  return(the_matrix)
})

# remove all genes with almost 0 counts
# threshold <- ncol(x)*2
mat_filt <- lapply(mat_raw, function(x) x[rowSums(x) > ncol(x)*2,])
mat_log <- lapply(mat_filt, log2)

# change -Inf values to 0
mat_log <- lapply(mat_log, function(x){x[x==-Inf] <- 0; return(x)})

# normalise all samples within a group to the sample with the largest read count
mat_norm <- lapply(mat_log, function(x){
  totalCounts <- colSums(x)
  corrections <- max(totalCounts)/totalCounts
  correctedCounts <- sweep(x, MARGIN=2, corrections, '*')
})


#=================================================
# calculates z-scores by using local distribution
#=================================================
zScores <- function(values_1,values_2, deviation_method="standard", slice_size=500) {
  
  average_values <- (values_1+values_2)/2
  
  order(average_values) -> sorted_indices
  
  order(sorted_indices) -> reverse_lookup
  
  sapply(1:length(values_1), function(x) {
    
    # if the 2 values are the same, there is no point calculating the z-score
    if((values_1[x] - values_2[x] == 0)){       
      z <- 0
      return(z)     
    }
    
    else{
      start <- reverse_lookup[x]-(slice_size/2)
      if (start < 0) start <- 0
      end <- start+slice_size
      if (end > length(values_1)) {
        end <- as.numeric(length(values_1))
        start <- end-slice_size
      }
      
      local_diffs <- as.double(values_1[sorted_indices[start:end]]-values_2[sorted_indices[start:end]])
      
      if(deviation_method=="standard"){
        
        # We assume a mean of 0 and calculate the sd
        local_dev <- sqrt(mean(local_diffs*local_diffs))
      }
      else if(deviation_method=="mad"){
        
        # median absolute deviation so that the standard deviation doesn't get totally skewed
        local_dev <- mad(local_diffs, center=0)
      }
      # again assuming a mean of 0
      z <- (values_1[x]-values_2[x]) / local_dev
    }
    return (z)
  })
}


# get some z-scores
z_scores <- lapply(mat_norm, function(x){
  
  if(ncol(x) != 2) stop("Each group must contain 2 samples")
  zScores(x[,1], x[,2])
})


# add in the z-scores, mean values and gene names to data frame so we've got all the information
annotated_df <- mapply(mat_norm, df, SIMPLIFY = FALSE, FUN=function(mat, df){
  
  if(ncol(mat) != 2) stop("Each group must contain 2 samples")
  #browser()
  gene_id <- rownames(mat)
  gene_name <- df$gene_name[match(gene_id,df$ensembl_id)]
  mean <- rowMeans(mat)
  z_score <- zScores(mat[,1], mat[,2])
  
  df <- data.frame(gene_id, gene_name, mat[,1], mat[,2], mean, z_score, row.names=NULL)
  
  #set the column names 
  colnames(df)[3:4] <- c(colnames(mat)[1], colnames(mat)[2])

  return(df)
})

# order data frame by absolute z-score so we can select the top 20
df_ordered <- lapply(annotated_df, function(x){
  
  x[order(abs(x$z_score), decreasing = TRUE),]
})


#----------------------
# text summary of data
#----------------------

# no of samples in sample sheet
text_file <- paste("../", names(samples)[1], "_summary.txt", sep="")
cat(paste(length(files), "files passed in\n"), file=text_file)
cat(paste(length(sample_sheet$sample), "sample names read from sample sheet\n"), file=text_file, append=TRUE)
cat(paste(length(matched_files), "sample names matched to files\n"), file=text_file, append=TRUE)

# no of raw genes
cat("\n\nnumber of raw genes\n", file=text_file, append=TRUE)
write.table(sapply(mat_raw, nrow), file = text_file, append = TRUE, col.names = FALSE, quote=FALSE)

# no of filtered genes
cat("\n\nnumber of filtered genes\n", file=text_file, append=TRUE)
write.table(sapply(mat_log, nrow), file = text_file, append = TRUE, col.names = FALSE, quote=FALSE)

#-------------------------------
# create pdf file for the plots
#-------------------------------
pdf(paste("../", names(samples)[1], ".pdf"))

par(mfrow= c(3,2))

# beanplot raw and normalised data
mapply(mat_log, mat_norm, names(mat_log), FUN=function(x,y,z){
  beanplot(as.data.frame(x), what=c(1,1,1,0), col="#6699ff", main=paste(z, "raw"))
  beanplot(as.data.frame(y), what=c(1,1,1,0), col="#6699ff", main=paste(z, "normalised"))
})

# z-score plot
mapply(df_ordered, names(df_ordered), FUN=function(x,y){
  plot(x$mean,x$z_score, pch=16, cex=0.8, col=densCols(x$mean,x$z_score, colramp = colorRampPalette(c("blue2","green","red"))), main = y)
  
  plot(x$mean,x$z_score, pch=16, cex=0.8, col="grey", main = y)
  points(x$mean[1:200],x$z_score[1:200], pch=16, cex=0.8, col="red")
})

dev.off()

#---------------------------------

# select top 200 genes by z-score to write out
# these can then be run through gprofiler
filtered <- lapply(df_ordered, head, n=200)

file_names_out <- paste(names(filtered), "_genes_top200_zscores.txt", sep="")

mapply(filtered, file_names_out, FUN=function(x,y) write.table(x$gene_name, file=y, row.names=FALSE,col.names = FALSE, quote=FALSE))

gene_names <- sapply(filtered, `[[`, "gene_name")
gene_names <- table(unlist(gene_names))
gene_names <- gene_names[gene_names>=2]
gene_names <- gene_names[order(gene_names, decreasing = TRUE)]

# write out genes that appear in >20% of the datasets
no_of_groups <- length(filtered)
threshold <- ceiling(no_of_groups*0.2)

gene_names_filt <- gene_names[gene_names>threshold]

message(paste(length(gene_names_filt), "genes appeared >", threshold, "times"))

write.table(gene_names_filt, file="most_frequently_occurring_genes.txt", quote=FALSE, col.names = FALSE, row.names=FALSE, sep="\t")
