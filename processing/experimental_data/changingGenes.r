# pass in a file containing a list of genes and associated values for 2 samples.
# Find the genes that change the most.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  
  filename <- args[1]
  output_filename <- gsub(pattern = ".txt", "_z_scores.txt", filename)
} else if (length(args)==2) {
  
  filename <- args[1]
  output_filename <- args[2]
}

file <- read.delim(filename)

# check all the gene ids are identical
sum(file$gene_id != file$gene_id.1)
    
counts <- file[,grepl("SRR", colnames(file))]
rownames(counts) <- file$gene_id

# remove any rows with 0 counts in both samples
counts  <- counts[rowSums(counts)!=0,]

# read count correction
correctForTotalReadCount <- function(df.counts){
  
  totalCounts <- colSums(df.counts)
  corrections <- max(totalCounts)/totalCounts
  sweep(df.counts, MARGIN=2, corrections, '*')
}


intensity.difference <- function (values.1,values.2) {
  average.values <- (values.1+values.2)/2
  
  order(average.values) -> sorted.indices
  
  order(sorted.indices) -> reverse.lookup
  
  sapply(1:length(values.1), function(x) {
    
    # we have a problem when we have a load of 0 values as when we have a sd of 0, the p value is 1, then when we do 1-local.p 
    #it's converted to 0. If all are set to 0.5, it messes up the q values, so I'm setting them to 1, it seems to work ok.....
    if((values.1[x] - values.2[x] == 0)){       
      local.p <- 1
      z <- 0
      return(c(z,local.p))      
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
      
      # assuming a mean of 0
      z <- ((values.1[x]-values.2[x]) - 0)/local.sd
      
      # Now we work out the p.value for the value we're actually looking at in the context of this distibution    
      pnorm(values.1[x]-values.2[x],mean=0,sd=local.sd) -> local.p
      
      if (local.p > 0.5){
        local.p <- (1 - local.p)
      }
      res <- c(z,local.p)
      
    }    
    return (res)
  }
  )
}

# log transform the data as we get some huge read counts which messes up the z-scores
counts[counts < 1] <- 1 # we don't want -Inf values 
counts <- log2(counts)

corrected_counts <- correctForTotalReadCount(counts)

x <- intensity.difference(corrected_counts[,1], corrected_counts[,2])
corrected_counts_stats <- cbind(corrected_counts, t(x))
colnames(corrected_counts_stats)[3:4] <- c("z","p")
corrected_counts_stats$abs_z <- abs(corrected_counts_stats$z)

# get the top 200 different genes
ordered_z <- corrected_counts_stats[order(corrected_counts_stats$abs_z, decreasing = TRUE),]

write.table(ordered_z , file=output_filename, quote=FALSE, sep="\t")

out_plotname <- gsub("_z_scores.txt", "_scatter.png", output_filename)

png(out_plotname, type="cairo", res=90, pointsize=16)

plot(ordered_z[,1], ordered_z[,2], pch=16, col="grey60", cex=0.7, xlab=colnames(ordered_z[1]), ylab =colnames(ordered_z[2]))
points(ordered_z[1:200,1], ordered_z[1:200,2], col="red", pch=16, cex=0.7)

dev.off()





