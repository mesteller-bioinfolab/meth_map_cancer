#!/usr/bin/env Rscript

# dependencies
library(multicore)

# set wd
setwd("/home/user/bsdata/hmr/")

# get data
dat <- read.delim("normal.all.hmr.counts.bed", head=F, stringsAsFactors=F)
names(dat) <- c("chr", "start", "end", "n")
dat$chr <- factor(dat$chr, levels=c(1:22, "X"))

# naive segmentation (select regions with more than min.samples and more than min.cpg)
hmr.pos <- function(dat, chr, min.samples=3, min.cpg=3){
      # subset of the data
      n <- dat$n[dat$chr == chr]
          names(n) <- dat$end[dat$chr == chr]

          # check where are the evicence of HMR
          hmr <- n >= min.samples
          #create codes
          codes <- names(n)
          names(codes) <- 1:length(n)
          #compute run lenghts
          runs <- rle(hmr)
          #set indexes as list per run
          ii <- tapply(1:sum(runs$lengths), rep(1:length(runs$lengths), runs$lengths), I)
          # select runs by size and status
          selected <- with(runs, which(lengths >= min.cpg & values))

          # translate selection into indexes
          out <- ii[selected]

          # get limits of the hmr
          out <- do.call(rbind, lapply(out, range))

          # convert indexes into positions
          start <- as.numeric(codes[as.character(out[, 1])]) - 1
          end <- as.numeric(codes[as.character(out[, 2])])

          # return bed format
          data.frame(chr=chr, start=start, end=end)
    }

# do the monkey!
res <- mclapply(unique(dat$chr), function(x) hmr.pos(dat, x, 7, 3))
res <- do.call(rbind, res)

# save output
write.table(res, "normal.common.hmr.bed", row.names=F, col.names=F, quote=F, sep="\t")

# do the monkey!
res <- mclapply(unique(dat$chr), function(x) hmr.pos(dat, x, 3, 3))
res <- do.call(rbind, res)

# save output
write.table(res, "normal.frequent.hmr.bed", row.names=F, col.names=F, quote=F, sep="\t")
