#!/usr/bin/env Rscript

# get arguments
input <- commandArgs(trailingOnly = TRUE)
n <- length(input)
smp1 <- input[1]
smp2 <- input[2]
outFile <- input[3]
chroms <- input[4:n]

# Usage message
if( n < 4){
  cat("\n")
  cat("Takes two tabix smoothed methylation files, an output dmr file and a set of chromosomes\n")
  cat("  ... and returns a dmr file\n")
  cat("\n")
  cat("    dmr file: (chr start end nCpG sign)\n")
  cat("\n")
  cat("    Usage: ./smooth.dmr.r sample1.smoothed.bed.gz sample2.smoothed.bed.gz outfile.bed chr1 chr2 chr3 ...\n")
  cat("\n")
  stop("Not enough arguments")
}


# dependencies
library(parallel)

# function to read data into BSseq objects
read.CGmeth <- function(methfile, chrom){
  # set temp file name
  tmp <- paste("/tmp/", paste(sample(c(0:9, letters), 12, rep=T), collapse=""), sep="")
  
  # extract selected chromosome
  system(paste("tabix", methfile, chrom, ">", tmp, sep=" "))
  
  # read data
  dat <- read.delim(tmp, head=F, stringsAsFactors=F, na.string=".")[, c(2, 5:9)]
  names(dat) <- c("pos", "u", "m", "s", "sl", "su")
  
  # remove tmp file
  system(paste("rm", tmp))
  
  # return
  dat
}

# funciton to look for dmr
naive.dmr <- function(x, pos){
  # compute run lenghts 
  runs <- rle(x)
  # set indexes as list per run
  ii <- tapply(1:sum(runs$lengths), rep(1:length(runs$lengths), runs$lengths), I)
  # select runs by size and status
  selected <- with(runs, which(lengths > 5 & values))
  # translate selection into indexes
  out <- as.data.frame(t(sapply(ii[selected], range)))
  # translate selection into positions
  out[,1] <- pos[out[,1]]
  out[,2] <- pos[out[,2]]
  out[,3] <- sapply(ii[selected], length)
  # return
  out
}

wrap <- function(smp1, smp2, chrom){
  # read data
  dat1 <- read.CGmeth(smp1, chrom)
  dat2 <- read.CGmeth(smp2, chrom)
  
  # search DMR
  hypo <- cbind(chr=chrom, naive.dmr(dat1$su < dat2$sl, dat1$pos), sense=-1)
  hyper <- cbind(chr=chrom, naive.dmr(dat1$sl > dat2$su, dat1$pos), sense=1)

  # merge and sort
  out <- rbind(hypo, hyper)
  out <- out[order(out[,1]),]
  gc()
  
  # return
  out
}


# run!
res <- mclapply(chroms, function(x) wrap(smp1, smp2, x), mc.cores=length(chroms))

# arrange and save results
res <- do.call(rbind, res)

write.table(res, outFile, row.names=F, col.names=F, quote=F, sep="\t")

