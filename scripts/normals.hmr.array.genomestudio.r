# dependencies
library(multicore)

# get array data
load("/home/user/bsdata/array/normals.genomestudio.RData")

############################
####    Tissue specific HMRs
############################

# get HMRs info
hmr <- read.delim("/home/user/bsdata/hmr/all.hmr.cgs.tab", head=F, stringsAsFactors=F)
names(hmr) <- c("cg", "hmr")

# get HMRs info
all.hmr <- read.delim("/home/user/bsdata/hmr/all.hmr.bed", head=F, stringsAsFactors=F)[,4]

# subset array data with probes in HMRs
dim(res)
res <- res[hmr$cg,]
dim(res)

# add hmr info
dim(res)
res <- merge(res, hmr)
dim(res)

# aggregate per hmr
resagg <- aggregate(res[, -c(1, ncol(res))], list(hmr=res$hmr), mean, na.rm=T)

# aggregate per hmr, requiere at least 5 cg per HMR
resagg.strict <- aggregate(res[, -c(1, ncol(res))], list(hmr=res$hmr),
                           function(x) {
                             out <- NA
                             if(length(x) > 4) out <- mean(x, na.rm=T)
                             out
                           })
resagg.strict <- resagg.strict[!apply(is.na(resagg.strict[, -1]), 1, all),]

# arrange and save objects
array.hmr.loose <- resagg[, -1]
rownames(array.hmr.loose) <- resagg$hmr
save(array.hmr.loose, file="/home/user/bsdata/array/array.hmr.loose.genomestudio.RData")

array.hmr.strict <- resagg.strict[, -1]
rownames(array.hmr.strict) <- resagg.strict$hmr
save(array.hmr.strict, file="/home/user/bsdata/array/array.hmr.strict.genomestudio.RData")

save(res, file="/home/user/bsdata/array/array.hmr.allcgs.genomestudio.RData")

# proportion of HMr with info at the array
numbers <- sapply(list(all=unique(all.hmr),
                       loose=rownames(array.hmr.loose),
                       strict=rownames(array.hmr.strict)),
                  function(x) sort(table(sapply(strsplit(x, "_"), "[", 1))))


sink("/home/user/bsdata/array/hmr.with.probes.genomestudio.txt")
numbers
sink()




