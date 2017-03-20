#!/usr/bin/env Rscript

# set wd
setwd("/home/user/bsdata/hmr")

# get sample info
info <- read.table("../samples.info.txt", stringsAsFactors=F)

normals <- info[substr(info$V3, 1, 4) == "Norm" & info$V2 != "Placenta" & info$V1 != "NB" & info$V1 != "103" & info$V1 != "Julia" & info$V1 != "W145", 1]

for(i in 1:length(normals)){
  ref <- normals[i]
  alt <- normals[-i]
  
  system(paste(paste("cat /home/user/bsdata/hmr/CpG_context_", ref, ".hmr.metrics.bed", sep=""), paste(" | bedtools intersect -v -a - -b /home/user/bsdata/hmr/CpG_context_", alt, ".hmr.bed ", sep="", collapse=""), " > /home/user/bsdata/hmr/", ref, ".hmr.specific.bed", sep=""))
}
