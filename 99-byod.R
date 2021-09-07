## https://www.embopress.org/doi/full/10.15252/msb.20188503
## Table EV1 from SI

library(readxl)
library(QFeatures)
library(tidyverse)

x <- read_xlsx("Table_EV1.xlsx", 4)


se <- readSummarizedExperiment(x, ecol = 3:31)
rownames(se) <- rowData(se)$Gene.ID

colnames(se)
longFormat(se[1:10, 1:5])

qf <- QFeatures(list(genes = se))
qf


filterFeatures(qf, ~ Gene.name %in% c("AFP", "NCOA6", "SLC8A3", "PLBD1", "C6orf47")) %>%
    longFormat() %>%
    as_tibble()

#####

data(hlpsms)
hl <- readQFeatures(hlpsms, ecol = 1:10, name = "psms")
hl$gr <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)

longFormat(hl[, , "proteins"]) %>%
    as_tibble() %>%
    full_join(data.frame(colData(hl)) %>% rownames_to_column("primary")) %>%
    group_by(gr, rowname) %>%
    summarise(mean = mean(value))
