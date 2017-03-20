# dependencies
library(multicore)
library(minfi)
options(stringsAsFactors=F)

# get cancer samples
cancer <- unique(read.delim("/home/user/bsdata/cancer_tissue_curated.txt", head=F, stringsAsFactors=F))
names(cancer) <- c("id", "tissue")

# get sebas sample sheet
sebas <- read.csv("/home/user/bsdata/array.samples.samplesheet.nopepe.csv")

# merge info
dat <- merge(cancer, sebas, by.x="id", by.y="Sample_Name", all.x=T)

# remove spaces in tissue names
dat$tissue <- gsub(" ", "_", dat$tissue)
# remove parenthesis in tissue names
dat$tissue <- gsub("(", "", dat$tissue, fixed=T)
dat$tissue <- gsub(")", "", dat$tissue, fixed=T)

# set wd
setwd("/home/user/minfi")

# function to store sample per tissue
wrap <- function(plate) {
  # get array data
  load(paste(plate, ".RData", sep=""))
  
  # set object name
  one <- getBeta(genomestudio)
  
  # remove original object
  #rm(genomestudio, swan)

  # get indexes of sampels in "dat"
  i <- match(paste(plate, colnames(one), sep="_"),
             paste(dat$Sample_Plate, dat$Sentrix_ID, dat$Sentrix_Position, sep="_"))

  # select samples present in dat
  one <- one[,which(!is.na(i)), drop=F]

  # set sample names
  colnames(one) <- na.omit(dat$id[i])

  # initialize output object
  out <- list()

  # for each tissue present in the plate
  for(ts in unique(na.omit(dat$tissue[i]))){
    # store all samples
    out[[ts]] <- one[,na.omit(dat$tissue[i]) == ts, drop=F]
  }

  # remove object
  rm(one)

  #output
  out
}


# store sample per tissue
res <- mclapply(unique(dat$Sample_Plate), wrap, mc.cores=10)
for(i in which(sapply(res, class) == "try-error")) res[[i]] <- wrap(unique(dat$Sample_Plate)[i])

# add mapping info
map <- read.delim("/home/user/resources/human/450k.manifest.1to22andX.onlypos.bed", head=F, stringsAsFactors=F)
names(map) <- c("chr", "start", "end", "cg")


# arrange data per tissue
out <- mclapply(unique(dat$tissue), function(ts) do.call(cbind, lapply(res, function(x) as.data.frame(x[[ts]])[map$cg,,drop=F])),
                mc.cores=10)
names(out) <- unique(dat$tissue)


# per tissue
names(map)[1] <- paste("#", names(map)[1], sep="")
mclapply(names(out), function(i){
write.table(cbind(map, out[[i]]),
            paste("/home/user/bsdata/array/tissue/",
                  "450k.betas.cancer.", i, ".genomestudio.bed", sep=""),
            quote=F, row.names=F, sep="\t")
})

# put all in a data frame
res <- lapply(res, function(x) do.call(cbind, x))
res <- do.call(cbind, res)
res <- as.data.frame(res)
res$cg <- rownames(res)

# save as a RData
save(res, file="/home/user/bsdata/array/cancers.genomestudio.RData")
