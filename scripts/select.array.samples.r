# dependencies
library(gdata)
options(stringsAsFactors=F)

# get sample sheet
info <- read.xls("/home/data/Samples on Array.xlsx", skip=6)

# get selected normal and cancer samples
normals <- read.delim("/home/user/bsdata/normal_tissue.txt", head=F)
cancers <- read.delim("/home/user/bsdata/cancer_tissue.txt", head=F)
names(normals) <- names(cancers) <- c("Sample_Name", "Tissue")
normals$status <- "normal"
cancers$status <- "cancer"
selected <- rbind(normals, cancers)

# merge info
dat <- merge(info, selected, all.y=T)

##  Sample names without info

# get samples without info
noinfo <- dat$Sample_Name[is.na(dat$Sentrix_ID)]

# remove that from data
prueba <- dat[!(dat$Sample_Name %in% noinfo),]

# get selected
selected.noinfo <- selected[selected$Sample_Name %in% noinfo,]

# change "-" to "_" and viceversa
selected.noinfo$Sample_Name <- gsub("-", "@", selected.noinfo$Sample_Name)
selected.noinfo$Sample_Name <- gsub("_", "-", selected.noinfo$Sample_Name)
selected.noinfo$Sample_Name <- gsub("@", "_", selected.noinfo$Sample_Name)

# merge info
selected.noinfo <- merge(info, selected.noinfo, all.y=T)

# add to previous
prueba <- rbind(prueba, selected.noinfo)
dat <- prueba
dat$Sample_Name[is.na(dat$Sentrix_ID)]

##  Duplicated sample names

# get duplicated IDs
dupIDs <- names(table(dat$Sample_Name))[table(dat$Sample_Name) > 1]

# curation
prueba <- dat
for(i in dupIDs){
  keep.row <- 0
  dup.data <- prueba[prueba$Sample_Name == i,]
  while(!(keep.row %in% rownames(dup.data))){
    cat("\n")
    print(dup.data[, c("Sample_Name", "Investigator", "Type", "Tissue")])
    cat("\n")
    keep.row <- readline("Which rownumber do you want to keep? ")
    cat("\n")
  }
  prueba <- prueba[!(rownames(prueba) %in% setdiff(rownames(dup.data), keep.row)),]
}

dat <- prueba
names(table(dat$Sample_Name))[table(dat$Sample_Name) > 1]

##  Remove centenarian and newborn blood samples
dat <- dat[!(dat$Tissue %in% c("Blood (old)", "Blood (newborn)")),]

##  Remove centenarian CD19
dat <- dat[dat$Sample_Name != "CD19_103_YEARS",]


## Remove Pepe Roman Samples
dat.nopepe <- dat[grep("pepe", dat$Coment, ignore.case=T, invert=T),]


## Total numbers
numbers <- rbind(table(dat$Tissue, dat$status), Total=apply(table(dat$Tissue, dat$status), 2, sum))
numbers.nopepe <- rbind(table(dat.nopepe$Tissue, dat.nopepe$status),
                        Total=apply(table(dat.nopepe$Tissue, dat.nopepe$status), 2, sum))

# save annotation files
write.csv(dat, "/home/user/bsdata/array.samples.samplesheet.csv", row.names=F)

write.csv(dat.nopepe, "/home/user/bsdata/array.samples.samplesheet.nopepe.csv", row.names=F)

# save selected samples
cancers <- dat.nopepe[dat.nopepe$status == "cancer", c("Sample_Name", "Tissue"),]
cancers <- cancers[order(cancers$Tissue, cancers$Sample_Name),]

normals <- dat.nopepe[dat.nopepe$status == "normal", c("Sample_Name", "Tissue"),]
normals <- normals[order(normals$Tissue, normals$Sample_Name),]

write.table(cancers, "/home/user/bsdata/cancer_tissue_curated.txt",
            quote=F, col.names=F, row.names=F, sep="\t")
write.table(normals, "/home/user/bsdata/normal_tissue_curated.txt",
            quote=F, col.names=F, row.names=F, sep="\t")


