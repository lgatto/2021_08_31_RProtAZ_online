library("QFeatures")
library("magrittr")
library("tidyverse")
library("limma")

## Homework: create a QFeatures object with the cptac data file
## msdata::quant(full.names = TRUE)
basename(f <- msdata::quant(full.names = TRUE))
names(read.delim(f))

i <- grep("Intensity\\.", names(read.delim(f)))

cptac_se <- readSummarizedExperiment(f, ecol = i,
                                     sep = "\t")

rowData(cptac_se)[1:10, 1:5]

anyDuplicated(rowData(cptac_se)[["Sequence"]])

readQFeatures(f, ecol = i,
              sep = "\t",
              fnames = "Sequence",
              name = "peptides")
QFeatures(list(peptides = cptac_se))

cptac_se <- readSummarizedExperiment(f, ecol = i,
                                     sep = "\t",
                                     fnames = "Sequence")

colData(cptac_se)

## Ex: populate the colData() with the experimental design of the
## data. See colnames(cptac_se)

colnames(cptac_se) <- sub("I.+\\.", "", colnames(cptac_se))

colData(cptac_se)$condition <- sub("_[7-9]", "", colnames(cptac_se))
cptac_se$id <- sub("^.+_", "", colnames(cptac_se))


colData(cptac_se)

names(rowData(cptac_se))

keep_var <- c("Sequence",
              "Proteins",
              "Leading.razor.protein",
              "PEP",
              "Score",
              "Reverse",
              "Potential.contaminant")

rowData(cptac_se) <- rowData(cptac_se)[, keep_var]

cptac_se <- zeroIsNA(cptac_se)
assay(cptac_se)

nNA(cptac_se)

MSnbase:::naplot(as(cptac_se, "MSnSet"),
                 col = c("black","white"))

barplot(nNA(cptac_se)$nNAcols$nNA)
table(nNA(cptac_se)$nNArows$nNA)

cptac_se <- filterNA(cptac_se, pNA = 4/6)

cptac_se

names(rowData(cptac_se))

table(rowData(cptac_se)$Reverse)

table(rowData(cptac_se)$Potential.contaminant)


cptac <- QFeatures(list(peptides = cptac_se))

colData(cptac) <- colData(cptac_se)

## Ex: remove contaminants and reverse hits with filterFeatures

cptac %>%
    filterFeatures(~ Reverse != "+")
