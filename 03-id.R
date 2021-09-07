px <- PXDataset("PXD000001")
rw <- pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")
sp <- Spectra(rw)

library(ggplot2)
library(dplyr)
library(PSM)
library(msdata)

calculateFragments("THSQEEMQHMQR")

## identification data
idf <- msdata::ident(full.names = TRUE)
id <- readPSMs(idf)

as_tibble(id) %>%
    ggplot(aes(x = MS.GF.RawScore,
               colour = isDecoy)) +
    geom_density()

id_filtered <- filterPSMs(id)


which(sp$spectrumId == "controllerType=0 controllerNumber=1 scan=2949")

table(table(id_filtered$spectrumID))

which(table(id_filtered$spectrumID) == 4)

i <- which(id_filtered$spectrumID == "controllerType=0 controllerNumber=1 scan=5490")

DT::datatable(as.data.frame(id_filtered[i, ]))

id_filtered <- reducePSMs(id_filtered, id_filtered$spectrumID)

which(id_filtered$spectrumID == "controllerType=0 controllerNumber=1 scan=5490")

id_filtered[1903, "modLocation"]

## combine raw and identification data

sp <- joinSpectraData(sp,
                      id_filtered,
                      by.x = "spectrumId",
                      by.y = "spectrumID")


spectraVariables(sp)


all(is.na(filterMsLevel(sp, 1)$MS.GF.RawScore))

table(is.na(filterMsLevel(sp, 2)$MS.GF.RawScore))


i <- which(sp$MS.GF.RawScore > 100)[1]

plotSpectra(sp[i])

sp$sequence[i]

addFragments(sp[i])

plotSpectra(sp[i], labels = addFragments)

as_tibble(peaksData(sp[i])[[1]])


## Unsupervised M/Z delta QC
## https://gist.github.com/lgatto/7de0bc9fd712b01604ef67c714580e78
## also available in Spectra version >= 1.3.8, that can be installed
## with BiocManager::install("RforMassSpectrometry/Spectra")

d <- computeMzDeltas(sp)
plotMzDelta(d)

## Homework:
##
## 1. plot the score densities for the 3 mzID files of the PXD022816
##    dataset.

PXD022816 <- PXDataset("PXD022816")

## 3 first mzML files
mzmls <- head(grep("mzML", pxfiles(PXD022816), value = TRUE), n = 3)
mzmls <- pxget(PXD022816, mzmls)

## 3 first msID files
mzids <- head(grep("mzID", pxfiles(PXD022816), value = TRUE), n = 3)
mzids <- pxget(PXD022816, mzids)

id <- readPSMs(mzids)
id

as_tibble(id) %>%
    ggplot(aes(x = MetaMorpheus.score,
               colour = isDecoy)) +
    geom_density() +
    facet_wrap(~ sub("^.+QEP", "QEP", id$spectrumFile))


## 2. combine the raw and identification data for the 3 MS
##    acquisitions of PXD022816. Hint: make sure you match PSMs to
##    scans in the right file.

sp <- Spectra(mzmls)

## primary key for spectra
sp$pkey <-
    paste0(sub("^.+_QEP", "QEP", basename(dataOrigin(sp))),
           gsub("^.+=", "::", sp$spectrumId))
head(sp$pkey)

## primary key for PSMs
id$pkey <-
    paste0(gsub("^.+\\QEP", "QEP", id$spectrumFile),
           sub("^.+=", "::", id$spectrumID))
head(id$pkey)

id_filtered <- filterPSMs(id)
id_filtered

sp <- joinSpectraData(sp, DataFrame(id_filtered), by.x = "pkey")
