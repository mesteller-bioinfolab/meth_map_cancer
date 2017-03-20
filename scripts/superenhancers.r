#!/usr/bin/env Rscript

options(stringsAsFactors=F)

# get arguments
input <- commandArgs(trailingOnly = TRUE)
outfile <- input[1]
files.bases <- input[-1]

# read data
dat.bases <- lapply(files.bases, function(x){
      dat <- read.delim(x, head=F)
          names(dat) <- c("chr", "start", "end", "enhancer", "signal", "noise", "overlap")
          dat$size <- dat$end - dat$start
          dat$overlap <- dat$overlap/dat$size
          dat$sample <- strsplit(x, ".", fixed=T)[[1]][3]
          dat$thres <- strsplit(x, ".", fixed=T)[[1]][5]
          dat$id <- paste("SE", dat$chr, dat$start, sep="_")
          dat
    })
dat.bases <- do.call(rbind, dat.bases)

# set information about superenhancers
dat.bases.info <- unique(dat.bases[, c(1:6)])

# get methylation info
dat.bases.thres100 <- data.frame(unclass(xtabs(overlap ~ enhancer + sample, subset(dat.bases, thres == "thres100"))))
dat.bases.thres100$enhancer <- rownames(dat.bases.thres100)
dat.bases.thres100 <- merge(dat.bases.info, dat.bases.thres100)

# arrange columns
dat.bases.thres100$chr <- gsub("chr", "", dat.bases.thres100$chr)
dat.bases.thres100 <- dat.bases.thres100[, c(2:4, 1, 5:ncol(dat.bases.thres100))]
names(dat.bases.thres100)[1] <- paste("#", names(dat.bases.thres100)[1], sep="")

# save all of them
write.table(dat.bases.thres100,
            outfile,
            sep="\t", quote=F, row.names=F)

