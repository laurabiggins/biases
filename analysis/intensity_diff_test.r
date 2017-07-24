# Rscript /bi/group/bioinf/Laura_B/bias_analysis/TIDIED/scripts/intensity_diff_test.r  Encode_CSHL_sample_sheet.txt encode_rna_seq_analysis/*gene_names.txt

args <- commandArgs(trailingOnly = TRUE)

#print(args[1])

#print ("=================")

if (length(args)<=1) {
  stop("Sample sheet and files to be processed must be supplied", call.=FALSE)
} else if (length(args)>1) {
  
  sample.sheet <- args[1]
  
  files <- args[2:length(args)]  
}

# print(sample.sheet)

# read in the sample sheet
sample.sheet <- read.delim(sample.sheet)

# split by gsm number
srrs <- split(sample.sheet$srr, f = sample.sheet$gsm)
ordered.gsms <- unique(sample.sheet$gsm)
gsm.description <- unique(sample.sheet$description)
gsm.info <- data.frame(gsm = ordered.gsms, description = gsm.description)

print(gsm.info)


# extract the SRR number from the file name so we can match it to the sample sheet 
files.split <- strsplit(files, split="/", fixed=TRUE)
file.names <- sapply(files.split, tail, n=1)
regex.matches <- regexpr("(SRR([0-9]+)_)",file.names)

file.names.srr <- substr(file.names, regex.matches, attr(regex.matches,"match.length")-1)
# remove empty fields
# file.names.srr <- file.names.srr[file.names.srr!=""]

#print(file.names.srr)



# we're just going to use the datasets that are in the sample sheet so we'll check that they match the 
# files passed in and use only the ones that overlap

srr.vector <- as.character(unlist(srrs))

# file names passed in that are not found in the sample sheet
print(paste(sum(!file.names.srr %in% srr.vector)," file names passed in that were not found in the sample sheet"))
print(file.names.srr[!file.names.srr %in% srr.vector])

# srrs from the sample sheet that are not found in the file names passed in
print(paste(sum(!srr.vector %in% file.names.srr)," srrs from the sample sheet that were not found in the file names passed in"))
print(srr.vector[!srr.vector %in% file.names.srr])

files.to.import <- files[file.names.srr %in% srr.vector]

# import datasets
datasets <- lapply(files.to.import, read.delim)



#############
# functions 
#############

# collapse list of datasets (which should be replicates) into a dataframe
collapseToDF <- function(datasets, dataset.names){
 
  # combine the files into 1 data frame
  sum.of.mismatches <- 0  
  
  # The ensembl ids should be the same for all the files
  for (i in 1:length(datasets))  for(j in 2:length(datasets)){
    sum.of.mismatches <- sum.of.mismatches + sum(datasets[[i]]$ensembl_id != datasets[[j]]$ensembl_id)    
  }
  if(sum.of.mismatches > 0){
    stop("ensembl ids don't match")
  }
  else{
    df.all <- data.frame(ensembl.id = datasets[[1]]$ensembl_id, gene.name = datasets[[1]]$gene_name, sapply(datasets, `[[`, 'score'))
    
    # remove the genes that have 0 counts
    df <- df.all[rowSums(df.all[,3:ncol(df.all)]) > 0, ]
    
  }  
  names(df)<- c("ensembl.id", "gene.name", as.character(dataset.names))
  
  return(df)
} 

getSignificantGenes <- function(dataset.values, gene.names, ensembl.ids, q.cutoff=0.05){
  # get the p values
  p.values <- intensityDiffMultipleSamples(dataset.values, ids = ensembl.ids)
  
  # do multiple testing correction
  q.values <- sapply(p.values, function(x){p.adjust(as.numeric(x), method = "BH")})
  
  # just select the rows/genes with q values < 0.05
  rows.i <- (rowSums(q.values < q.cutoff)) > 0
  
  selected <- q.values[rows.i, ]
  
  genes <- gene.names[rows.i]
  
  return(data.frame(ids = genes, selected))
  
}

intensityDiffMultipleSamples <- function(df, ids){
  
  results <- data.frame(row.names = ids)
  columns.to.drop <- c("ensembl.id", "gene.name")
  df.scores <- df[, !(names(df) %in% columns.to.drop)]
  
  for (i in 1:(ncol(df.scores)-1))  for(j in (i+1):ncol(df.scores))  if(i != j){
    
    textString <- paste("now comparing ", names(df.scores)[i], " to ", names(df.scores)[j], sep = "")
    print(textString)
    
    p.vals <- intensity.difference(df.scores[,i], df.scores[,j])
    comparison.name <- paste(colnames(df.scores)[i], colnames(df.scores)[j], sep = "_")
    results[,comparison.name] <- p.vals
    
  }  
  return(results)
}


# This is the intensity difference function that is used in SeqMonk
intensity.difference <- function (values.1,values.2) {
  
  # log transform
  v1 <- log2(values.1)
  v2 <- log2(values.2)
  
  v1[v1==-Inf] <- 0
  v2[v2==-Inf] <- 0
  
  values.1 <- v1
  values.2 <- v2
  
  average.values <- (values.1+values.2)/2
  
  order(average.values) -> sorted.indices
  
  order(sorted.indices) -> reverse.lookup
  
  sapply(1:length(values.1), function(x) {
    
    # we have a problem when we have a load of 0 values as when we have a sd of 0, the p value is 1, then when we do 1-local.p 
    #it's converted to 0. If all are set to 0.5, it messes up the q values, so I'm setting them to 1, it seems to work ok.....
    if((values.1[x] - values.2[x] == 0)){       
      local.p <- 1
      return(local.p)     
    }
    
    else{
      
      start <- reverse.lookup[x]-250
      if (start < 0) start <- 0
      end <- start+500
      if (end > length(values.1)) {
        end <- as.numeric(length(values.1))
        start <- end-500
      }
      
      local.diffs <- as.double(values.1[sorted.indices[start:end]]-values.2[sorted.indices[start:end]])
      
      # We assume a mean of 0 and calculate the sd
      sqrt(mean(local.diffs*local.diffs)) -> local.sd
      
      # Now we work out the p.value for the value we're actually looking at in the context of this distibution    
      pnorm(values.1[x]-values.2[x],mean=0,sd=local.sd) -> local.p
      
      if (local.p > 0.5){
        local.p <- (1 - local.p)
      } 
    }    
    return (local.p)
  }
  )
}


# this is not clever - if there is a really high no somewhere that will affect the correction
# just pass in the counts
correctForTotalReadCount <- function(df.counts){
  
  totalCounts <- colSums(df.counts)
  barplot(totalCounts, cex.names=0.7, las =2)
  
  maxCount <- max(totalCounts)
  
  corrections <- maxCount/totalCounts
  
  correctedCounts <- sweep(df.counts, MARGIN=2, corrections, '*')
  
  barplot(colSums(correctedCounts), cex.names=0.7, las =2)
  
  return(correctedCounts)
}


# split by gsm number
srrs <- split(sample.sheet$srr, f = sample.sheet$gsm)
ordered.gsms <- unique(sample.sheet$gsm)
gsm.description <- unique(sample.sheet$description)
gsm.info <- data.frame(gsm = ordered.gsms, description = gsm.description)

print(gsm.info)


# extract the SRR number from the file name so we can match it to the sample sheet 
files.split <- strsplit(files, split="/", fixed=TRUE)
file.names <- sapply(files.split, tail, n=1)
regex.matches <- regexpr("(SRR([0-9]+)_)",file.names)

file.names.srr <- substr(file.names, regex.matches, attr(regex.matches,"match.length")-1)
# remove empty fields
# file.names.srr <- file.names.srr[file.names.srr!=""]

#print(file.names.srr)


# we're just going to use the datasets that are in the sample sheet so we'll check that they match the 
# files passed in and use only the ones that overlap

srr.vector <- as.character(unlist(srrs))

# file names passed in that are not found in the sample sheet
print(paste(sum(!file.names.srr %in% srr.vector)," file names passed in that were not found in the sample sheet"))
print(file.names.srr[!file.names.srr %in% srr.vector])

# srrs from the sample sheet that are not found in the file names passed in
print(paste(sum(!srr.vector %in% file.names.srr)," srrs from the sample sheet that were not found in the file names passed in"))
print(srr.vector[!srr.vector %in% file.names.srr])


#=====================
# Importing the files
#=====================

# this uses lapply so we can deal with multiple gsms. It should also work with a single gsm.

files.to.import <- lapply(srrs, function(x) files[file.names.srr %in% as.character(x)])

# import datasets
datasets <- lapply(files.to.import, function(x) lapply(x, read.delim))

# the shortened names
dataset.names <- lapply(srrs, function(x) x[x %in% file.names.srr])

# collapse list of datasets (which should be replicates) into a dataframe
df <- mapply(collapseToDF, datasets, dataset.names, SIMPLIFY = FALSE)

# read count correction
countData <- lapply(df, function(x) correctForTotalReadCount(x[,3:ncol(x)]))    

# or no correction
#countData <- lapply(df, function(x) x[,3:ncol(x)])


res <- mapply(countData, df, FUN=function(x,y) getSignificantGenes(x, y[,"gene.name"],y[,"ensembl.id"]), SIMPLIFY = FALSE)

res.ordered <- lapply(res, function(x){
  #browser()
  # We'll order by the highest number of q-values that are < 0.05
  res.matrix <- as.matrix(x[,2:ncol(x)])
  rownames(res.matrix) <- x$ids
  
  freq <- apply(res.matrix, MARGIN=1, FUN = function(x)sum(x<0.05))
  x[order(freq, decreasing = TRUE),]
})


#===================
# write out results
#===================

# write out the list of genes and q values
mapply(res.ordered, names(res.ordered), FUN = function(x,y) write.table(x = x, file = paste(y,"_results_table.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE))

# just a list of genes
mapply(res.ordered, names(res.ordered), FUN = function(x,y) write.table(x = x$ids, file = paste(y,"_genelist.txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names=FALSE))






