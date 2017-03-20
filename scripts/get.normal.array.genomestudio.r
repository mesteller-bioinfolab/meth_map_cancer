# dependencies
library(multicore)
library(minfi)
options(stringsAsFactors=F)

# get normal samples
normal <- read.delim("/home/user/bsdata/normal_tissue_curated.txt", head=F)
names(normal) <- c("id", "tissue")

# get sebas sample sheet
sebas <- read.csv("/home/user/bsdata/array.samples.samplesheet.nopepe.csv")

# merge info
dat <- merge(normal, sebas, by.x="id", by.y="Sample_Name", all.x=T)

# remove spaces in tissue names
dat$tissue <- gsub(" ", "_", dat$tissue)
# remove parenthesis in tissue names
dat$tissue <- gsub("(", "", dat$tissue, fixed=T)
dat$tissue <- gsub(")", "", dat$tissue, fixed=T)

# set wd
setwd("/home/user/minfi")

# get selected samples
res <- mclapply(unique(dat$Sample_Plate), function(smet){

  # load batch
  load(paste(smet, ".RData", sep=""))
  betas <- getBeta(genomestudio)
  rm(genomestudio, swan)
  
  # get samples id 
  ids <- paste(dat$Sentrix_ID[dat$Sample_Plate == smet], dat$Sentrix_Position[dat$Sample_Plate == smet], sep="_")

  # set codes
  id.codes <- dat$id[dat$Sample_Plate == smet]
  names(id.codes) <- ids

  # detect absent samples
  NA.samples <- id.codes[setdiff(ids, colnames(betas))]

  # set present samples
  ids <- intersect(ids, colnames(betas))

  # subset the samples of interest
  out <- betas[, ids]

  # change sample names
  colnames(out) <- id.codes[colnames(out)]

  # set missing samples
  attr(out, "missing") <- NA.samples

  # output
  out
})

# backup
res.old <- res
              
# put all in a data frame
res <- do.call(cbind, res)
res <- as.data.frame(res)
res$cg <- rownames(res)

# save as a RData
save(res, file="/home/user/bsdata/array/normals.genomestudio.RData")

# add mapping info
map <- read.delim("/home/user/resources/human/450k.manifest.1to22andX.onlypos.bed", head=F, stringsAsFactors=F)
names(map) <- c("chr", "start", "end", "cg")

res <- merge(map, res)
res <- res[, c(2:4, 1, 5:ncol(res))]

res <- res[order(res$chr, res$start),]

# save results
# per file
mclapply(5:ncol(res), function(i){
  write.table(res[, c(1:4, i)],
              paste("/home/user/bsdata/array/",
                    "450k.betas.", names(res)[i], ".genomestudio.bed", sep=""),
              quote=F, row.names=F, col.names=F, sep="\t")
})

# per tissue
mclapply(unique(dat$tissue), function(i){
write.table(res[, c(1:4, na.omit(match(dat$id[dat$tissue == i], colnames(res))))],
            paste("/home/user/bsdata/array/tissue/",
                  "450k.betas.", i, ".genomestudio.bed", sep=""),
            quote=F, row.names=F, sep="\t")
})


