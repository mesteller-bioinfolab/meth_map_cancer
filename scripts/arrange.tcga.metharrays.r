# dependencies
library(multicore)
options(stringsAsFactors=F)

# set wd
setwd("/home/user/TCGA_validation")

# get CpG cpositions
map <- read.delim("/home/user/resources/human/450k.manifest.1to22andX.onlypos.bed", head=F, stringsAsFactors=F)
names(map) <- c("chr", "start", "end", "cg")
map <- map[substr(map$cg, 1, 2) == "cg",]

# get sample sheet
samplesheet <- read.csv("/home/data/Samplesheet_TCGA_allsamples.csv")
samplesheet$key <- paste(samplesheet$Sentrix_ID, samplesheet$Sentrix_Position, sep="_")

# get files
files <- dir(patt=".RData")

# arrange data as bed files
mclapply(files, function(file){
  # load data
  load(file)
  
  # filter out based on detection p-value
  betas <- as.matrix(betas)
  pval <- as.matrix(pval)
  betas[pval > 0.05] <- NA

  # change sample names
  colnames(betas) <- samplesheet$Sample_Name[match(colnames(betas), samplesheet$key)]
  
  # add genomic position info
  out <- cbind(map, betas[map$cg,])
  
  # write output
  if(grepl("_norm", file)){
    write.table(out,
                paste("/home/user/bsdata/tcga_meth/", strsplit(file, "_")[[1]][1], ".bed", sep=""),
                row.names=F, sep="\t", quote=F)
  }else{
    write.table(out,
                paste("/home/user/bsdata/tcga_meth/",
                      tolower(strsplit(file, ".", fixed=T)[[1]][1]), "_complete.bed", sep=""),
                row.names=F, sep="\t", quote=F)
  }
    
})
