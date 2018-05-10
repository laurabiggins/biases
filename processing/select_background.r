setwd("C:/Users/bigginsl/Desktop")
library("data.table")

all_genfo <- fread("Mus_musculus.GRCm38.80_gene_info.txt")

genelist <- scan("O:/Training/FAGL/Gene List Course Data/Presenting_results/genelist1.txt", what="character")

# check how many genes are not found in the background set
match_loc <- match(genelist, all_genfo$gene_name)
print(paste(sum(!is.na(match_loc)), "genes found in background"))
print(paste(sum(is.na(match_loc)), "genes not found in background....."))
print(genelist[is.na(match_loc)]) 

# filter for genes that are in the background and only pass these to the function
query_genes_in_bg <- genelist[genelist%in%all_genfo$gene_name]

getProbabilities <- function(wholeBackground, query){
  predict_density = approxfun(density(query, n=length(wholeBackground))) #function that approximates dens.obs
  prob <- predict_density(wholeBackground)
  prob[is.na(prob)] <- 0
  prob
}

matchBackground <- function(whole_background_genfo, query_genes, gene_length=TRUE, GC_content=TRUE, sample_size=5000, norm=TRUE){
  #browser()
  if(gene_length!=TRUE & gene_length!=FALSE){
    stop(paste("gene_length must be true or false, not ", gene_length))
  }
  if(GC_content!=TRUE & GC_content!=FALSE){
    stop(paste("GC_content must be true or false, not ", GC_content))
  }

  match_loc <- match(query_genes, whole_background_genfo$gene_name)
  query_genfo <- whole_background_genfo[na.omit(match_loc),]
  
  if(gene_length==TRUE & GC_content==TRUE){
    length_prob <- getProbabilities(log2(whole_background_genfo$length), log2(query_genfo$length))
    length_norm <- (length_prob-min(length_prob))/(max(length_prob)-min(length_prob))
    
    gc_prob <- getProbabilities(whole_background_genfo$GC_content, query_genfo$GC_content)
    gc_norm <- (gc_prob-min(gc_prob))/(max(gc_prob)-min(gc_prob))
    
    
    ifelse(norm==TRUE, geo_mean <- sqrt(length_norm*gc_norm), geo_mean <- sqrt(length_prob*gc_prob))
    
    probabilities <- geo_mean
    
  }
  else if(gene_length==TRUE & GC_content==FALSE){
    length_prob <- getProbabilities(log2(whole_background_genfo$length), log2(query_genfo$length))
    probabilities <- length_prob
  }  
  else if(GC_content==TRUE & gene_length==FALSE){
    gc_prob <- getProbabilities(whole_background_genfo$GC_content, query_genfo$GC_content)
    probabilities <- gc_prob
  }
  
  genes <- sample(whole_background_genfo$gene_id, sample_size, replace=FALSE, prob=probabilities)
  whole_background_genfo[whole_background_genfo$gene_id %in% genes]
}


# It is difficult to get multiple distributions to match well, but this seems to be a reasonable method.
# This method does work well when the probabilities are 0 for a condition - those genes will not be selected as the geometric mean is used so 0 * x = 0
# 



