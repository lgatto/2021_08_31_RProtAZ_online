## install from CRAN
install.packages("BiocManager")

library("BiocManager")

## install from Bioc
BiocManager::install("msdata")
BiocManager::install("mzR")
BiocManager::install("rpx")
BiocManager::install("ProtGenerics")
BiocManager::install("MsCoreUtils")
BiocManager::install("QFeatures")
BiocManager::install("Spectra")

## install from github
BiocManager::install("RforMassSpectrometry/PSM")
BiocManager::install("lgatto/rpx")

## Bioc version
BiocManager::version()

## R version
base::version

## rpx package version
packageVersion("rpx")
