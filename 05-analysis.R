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
cptac <- cptac %>%
    filterFeatures(~ Reverse != "+") %>%
    filterFeatures(~ Potential.contaminant != "+") %>%
    filterFeatures(~ PEP < 0.05)

limma::plotDensities(assay(cptac[[1]]))

cptac <- logTransform(cptac,
                      i = "peptides",
                      name = "log_peptides")


limma::plotDensities(assay(cptac[["log_peptides"]]))


## Ex: normalise the data
cptac <- normalize(cptac,
                   i = "log_peptides",
                   name = "lognorm_peptides",
                   method = "center.median")

plot(cptac)

BiocManager::install("limma")
BiocManager::install("MSnbase")

par(mfrow = c(1, 3))
limma::plotDensities(assay(cptac[["peptides"]]))
limma::plotDensities(assay(cptac[["log_peptides"]]))
limma::plotDensities(assay(cptac[["lognorm_peptides"]]))

## Ex: aggregate peptide features into proteins
cptac <- aggregateFeatures(cptac,
                           i = "lognorm_peptides",
                           name = "proteins_med",
                           fcol = "Leading.razor.protein",
                           fun = colMedians, ## try MsCoreUtils::robustSummary
                           na.rm = TRUE)


cptac <- aggregateFeatures(cptac,
                  i = "lognorm_peptides",
                  name = "proteins_robust",
                  fcol = "Leading.razor.protein",
                  fun = MsCoreUtils::robustSummary,
                  na.rm = TRUE)

BiocManager::install("RforMassSpectrometry/QFeatures")

table(rowData(cptac[["proteins_med"]])$.n)

library("factoextra")

pca_pep <-
    cptac[["lognorm_peptides"]] %>%
    filterNA() %>%
    assay() %>%
    t() %>%
    prcomp(scale = TRUE, center = TRUE) %>%
    fviz_pca_ind(habillage = cptac$condition)


pca_prot <-
    cptac[["proteins_med"]] %>%
    filterNA() %>%
    assay() %>%
    t() %>%
    prcomp(scale = TRUE, center = TRUE) %>%
    fviz_pca_ind(habillage = cptac$condition)

library(patchwork)
pca_pep + pca_prot



cptac["P02787ups|TRFE_HUMAN_UPS", ,
      c("lognorm_peptides", "proteins_med")] %>%
    longFormat() %>%
    as_tibble() %>%
    mutate(condition = ifelse(grepl("A", colname), "A", "B")) %>%
    ggplot(aes(x = colname,
               y = value,
               colour = rowname,
               shape = condition)) +
    geom_point(size = 3) +
    geom_line(aes(group = rowname)) +
    facet_grid(~ assay) +
    ggtitle("P02787ups|TRFE_HUMAN_UPS")
