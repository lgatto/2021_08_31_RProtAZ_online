## BiocManager::install("QFeatures")
library(QFeatures)

data(feat1)
feat1

feat1[[1]]
feat1[["psms"]]

colData(feat1)

colData(feat1)$X  <- c("X1", "X2")
feat1$Y  <- c("Y1", "Y2")
colData(feat1)

assay(feat1[["psms"]])
rowData(feat1[["psms"]])

## Feature aggregation

feat1 <- aggregateFeatures(feat1,
                           i = "psms",
                           fcol = "Sequence",
                           name = "peptides",
                           fun = colMeans)

feat1

rowData(feat1[["peptides"]])

assay(feat1[["peptides"]])

feat1 <- aggregateFeatures(feat1,
                           i = "peptides",
                           fcol = "Protein",
                           name = "proteins",
                           fun = colMedians)

assay(feat1[[3]])
rowData(feat1[[3]])


## Subsetting and filtering

pA <- feat1["ProtA", , ]
pV <- filterFeatures(feat1, ~ pval < 0.05)
## filter features that do not location to the MT
filterFeatures(feat1, ~ location != "Mitochondrion")


## How to create a QFeatures object
data(hlpsms)

hl <- readQFeatures(hlpsms, ecol = 1:10, name = "psms")
hl

hl <- aggregateFeatures(hl,
                        i = "psms",
                        fcol = "Sequence",
                        name = "peptides",
                        fun = colMeans)

hl <- aggregateFeatures(hl,
                        i = "peptides",
                        fcol = "ProteinGroupAccessions",
                        name = "proteins",
                        fun = colMeans)


hl

stat3 <- hl["P42227-2", , ]

stat3_df <- data.frame(longFormat(stat3))

stat3_df$assay <- factor(stat3_df$assay,
                         levels = c("psms", "peptides",
                                    "proteins"))

library(ggplot2)

ggplot(data = stat3_df,
       aes(x = colname,
           y = value,
           group = rowname)) +
    geom_line() + geom_point() +
    facet_grid(~ assay)


stat <- hl[c("P42227-2", "P42225"), , ]
stat

stat_df <- data.frame(longFormat(stat))
stat_df$stat3 <- ifelse(stat_df$rowname %in% stat3_df$rowname,
                        "STAT3", "STAT1")
stat_df$assay <- factor(stat_df$assay,
                        levels = c("psms", "peptides", "proteins"))

ggplot(data = stat_df,
       aes(x = colname,
           y = value,
           group = rowname)) +
    geom_line() + geom_point() +
    facet_grid(stat3 ~ assay)


## Homework: create a QFeatures object with the cptac data file
## msdata::quant(full.names = TRUE)
