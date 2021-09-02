## http://lgatto.github.io/rpx/
library("rpx")

px <- PXDataset("PXD000001")
px

## error
## PXDataset("PXDqhjkagdejs")

fls <- pxfiles(px)
fls

pxtax(px)

pxurl(px)

pxref(px)

fas <- pxget(px, fls[1])

pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")

rw <- pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")

## if you have your own files
rw <- "rawdata/myMSrun.mzML"


## advanced feature: delete project from cache
library(BiocFileCache)

## bfcremove(rpxCache(), "BFC1")
bfcremove(rpxCache(), "BFC1770")

library("msdata")

proteomics()
ident()
quant()

cptac_file <- quant(full.names = TRUE)
cptac_file


## conversion from binary to open format
## https://proteowizard.sourceforge.io/ (Windows only, also as a GUI)
## https://github.com/compomics/ThermoRawFileParser (all platforms)


## Homework: Download the 3 first mzML (raw data) and mzID
## (identification data) files from the PXD022816 project
## (Morgenstern, Barzilay, and Levin 2021).

px2 <- PXDataset("PXD022816")

## 3 first mzML files
mzmls <- head(grep("mzML", pxfiles(px2), value = TRUE), n = 3)
mzmls <- pxget(px2, mzmls)

## 3 first msID files
mzids <- head(grep("mzID", pxfiles(px2), value = TRUE), n = 3)
mzids <- pxget(px2, mzids)
