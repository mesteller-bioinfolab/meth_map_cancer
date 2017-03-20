# dependencies
library(ggplot2); theme_set(theme_bw())
library(scales)
library(multicore)
library(RColorBrewer)
options(stringsAsFactors=F)

# set wd
setwd("/home/user/bsdata/pca")

# get file names
files <- dir(patt="CpG.raw4cove")

# get data
dat <- mclapply(files, read.delim, head=F, na.string=".")

# get sample names
sample.names <- read.table("sample.names")[,1]

# remove OLC samples and CpG with missing
dat <- mclapply(dat, function(x){
  names(x) <- sample.names
  x <- x[, grep("OCL", names(x), invert=T)]
  na.omit(x)})

dat.old <- dat

# merge all chromosomes
dat <- do.call(rbind, dat)

# compute principal components centering the CpGs
pca <- princomp(t(scale(t(dat), scale=F)))

save(pca, file="../results/epigenomics.pca.RData")

# get loadings
ggdat <- as.data.frame(unclass(pca$loadings))

# get sample name
ggdat$sample <- names(dat)

# get sample info
info <- read.csv("../samples.info.newnames.csv", head=F)
names(info) <- c("sample", "tissue", "status", "type", "newname")
info$control <- "tumor"
info$control[grep("normal", info$status, ignore=T)] <- "normal"
ggdat <- merge(ggdat, info)

# proportion of variance explained
labels <- round(100*pca$sdev^2/sum(pca$sdev^2))
labels <- paste(names(labels), ": ", labels, "%", sep="")

# plots
pdf("../results/plots/epigenomics.pca.pdf")
ss <- 3
ggplot(ggdat,
       aes(x=Comp.1, y=Comp.2, label=newname, col=control)) +
  geom_text(size=ss, alpha=.85) +
  ggtitle("Control") +
  xlab(labels[1]) + ylab(labels[2])

ggplot(ggdat,
       aes(x=Comp.1, y=Comp.2, label=newname, col=type)) +
  geom_text(size=ss, alpha=.85) +
  ggtitle("Type") +
  xlab(labels[1]) + ylab(labels[2])

ggplot(ggdat,
       aes(x=Comp.1, y=Comp.2, label=newname, col=tissue)) +
  geom_text(size=ss, alpha=.85) +
  ggtitle("Tissue") +
  xlab(labels[1]) + ylab(labels[2])

ggplot(mapping=
       aes(x=Comp.1, y=Comp.2, label=newname, col=tissue)) +
  geom_point(data=subset(ggdat, newname %in% c("Prostate", "NB", "Y26", "Y103", "CD19")), size=ss, alpha=.7) +
  geom_text(data=subset(ggdat, ! (newname %in% c("Prostate", "NB", "Y26", "Y103", "CD19"))), size=ss) +
  ggtitle("Tissue") +
  xlab(labels[1]) + ylab(labels[2])

ggplot(ggdat,
       aes(x=Comp.1, y=Comp.2, label=newname, col=type)) +
  geom_text(size=ss, alpha=.85) +
  ggtitle("Type") +
  xlab(labels[1]) + ylab(labels[2]) +
  facet_wrap(~ tissue)

ggplot(ggdat,
       aes(x=Comp.1, y=Comp.2, label=newname, col=control)) +
  geom_text(size=ss, alpha=.85) +
  ggtitle("Control") +
  xlab(labels[1]) + ylab(labels[2]) +
  facet_wrap(~ tissue)


dev.off()

theplot <- ggplot(ggdat,
       aes(x=Comp.1, y=Comp.2)) +
  geom_point() +
  xlab(labels[1]) + ylab(labels[2]) +
  theme(axis.line=element_line(size=1))
ggsave("../results/plots/epigenomics.pca.points.pdf", theplot)

##  Clustering

names(dat) <- info$newname[match(names(dat), info$sample)]

tdat <- t(dat)

d.raw <- dist(tdat)
d.cent <- dist(scale(tdat, scale=F))
d.scale <- dist(scale(tdat))

clust.raw <- hclust(d.raw)
clust.cent <- hclust(d.cent)
clust.scale <- hclust(d.scale)

pdf("../results/plots/epigenomics.cluster.pdf")
plot(clust.raw, main="Raw")
plot(clust.cent, main="Center")
plot(clust.scale, main="Scale")
dev.off()


# Get SD of CpGs to filter out less variable ones
cpg.sd <- sqrt(rowSums((dat - rowMeans(dat))^2)/(ncol(dat) - 1))

summary(cpg.sd)

pdf("../results/plots/cpg.sd.pdf")
plot(density(cpg.sd))
abline(v=c(.2, .25), col="red")
hist(cpg.sd, 100)
abline(v=c(.2, .25), col="red")
dev.off()

# redo the clustering only with the ~ 50% most variant CpGs

tdat.filt <- t(dat[cpg.sd > .2,])

d.raw <- dist(tdat.filt)
d.cent <- dist(scale(tdat.filt, scale=F))
d.scale <- dist(scale(tdat.filt))

clust.raw <- hclust(d.raw)
clust.cent <- hclust(d.cent)
clust.scale <- hclust(d.scale)

pdf("../results/plots/epigenomics.cluster.filt.pdf")
plot(clust.raw, main="Raw")
plot(clust.cent, main="Center")
plot(clust.scale, main="Scale")
dev.off()


tdat.filt25 <- t(dat[cpg.sd > .25,])

d.raw <- dist(tdat.filt25)
d.cent <- dist(scale(tdat.filt25, scale=F))
d.scale <- dist(scale(tdat.filt25))

clust.raw <- hclust(d.raw)
clust.cent <- hclust(d.cent)
clust.scale <- hclust(d.scale)

pdf("../results/plots/epigenomics.cluster.filt25.pdf")
plot(clust.raw, main="Raw")
plot(clust.cent, main="Center")
plot(clust.scale, main="Scale")
dev.off()

