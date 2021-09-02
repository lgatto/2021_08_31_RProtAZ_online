## https://rformassspectrometry.github.io/Spectra
library("Spectra")
library("magrittr")

## a dummy example
spd <- DataFrame(
    msLevel = c(1L, 2L),
    rtime = c(10, 20))

spd$mz <- list(c(100, 220, 300),
               c(112, 134, 220, 300))
spd$intensity <- list(c(23, 345, 112),
                      c(72, 199, 234, 82))

sp0 <- Spectra(spd)

## accessing
spectraData(sp0)
spectraVariables(sp0)
peaksData(sp0)
peaksData(sp0)[[1]]
peaksData(sp0)[[2]]

## subsetting
sp0[1]
sp0[2]
sp0[c(1, 2, 2, 1, 1, 1)]

## a real example
px <- PXDataset("PXD000001")
pxfiles(px)
f <- pxget(px, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML")

sp <- Spectra(f)
sp

## - check the number of scan with length()
length(sp)
## - explore that data with the accessors above
spectraVariables(sp)
spectraData(sp[1:10])
peaksData(sp[1:10])

?Spectra

msLevel(sp)
sp$msLevel

length(msLevel(sp))

precursorMz(sp)
sp$precursorMz

centroided(sp)

filterMsLevel(sp, 2L) %>%
    precursorMz()

filterMsLevel(sp, 2L) %>%
    precursorIntensity()

## how many MS levels?
unique(msLevel(sp))

## are my data centroided or in profile mode?
unique(centroided(sp))

table(msLevel(sp), centroided(sp))

## What MS2 scan has the highest base peak intensity?
sp2 <- filterMsLevel(sp, 2L)
i <- which.max(sp2$basePeakIntensity)
sp2[i]

rtime(sp2[i])

## Plot it
plotSpectra(sp2[i])

## Filter on intensity
summary(intensity(sp2[i])[[1]])

mz(sp2[i])

plotSpectra(filterIntensity(sp2, 1522859)[i])

par(mfrow = c(2, 1))
plotSpectra(sp2[i])
plotSpectra(filterIntensity(sp2, 1522859)[i])

(fls <- dir(system.file("sciex",
                        package = "msdata"),
            full.names = TRUE))

sp_sciex <- Spectra(fls)

table(dataOrigin(sp_sciex))


## Peak picking (profile to centroided)

plotSpectra(sp[2807], xlim = c(521.2, 522.5),
            type = "l")

par(mfrow = c(2, 1))

plotSpectra(sp[2807], xlim = c(521.2, 522.5),
            type = "h")

pickPeaks(sp[2807]) %>%
    filterIntensity(1e7) %>%
    plotSpectra(xlim = c(521.2, 522.5))


## Visualisation exercices

## The chromatogram can be created by extracting the totIonCurrent and
## rtime variables for all MS1 spectra. Annotate the spectrum of
## interest (index 2807)

x <- filterMsLevel(sp, 1L) %>%
    spectraData(c("rtime", "totIonCurrent"))
plot(x$rtime, x$totIonCurrent,
     type = "l")
abline(v = rtime(sp[2807]), col = "red")

## The filterPrecursorScan() function can be used to retain parent
## (MS1) and children scans (MS2), as defined by an acquisition
## number. Use it to extract the MS1 scan of interest and all its MS2
## children.

ms_2 <- filterPrecursorScan(sp, 2807)

## Plot the MS1 spectrum of interest and highlight all the peaks that
## will be selected for MS2 analysis.

plotSpectra(ms_2[1], xlim = c(400, 1000))
abline(v = precursorMz(ms_2)[-1], col = "grey")

## Zoom in mz values 521.1 and 522.5 to reveal the isotopic envelope
## of that peak.

plotSpectra(ms_2[1], xlim = c(521.1, 522.5), type = "l")
abline(v = precursorMz(ms_2)[-1], col = "grey")

## The plotSpectra() function is used to plot all 10 MS2 spectra in
## one call.

plotSpectra(ms_2[-1])

## Homework:
##
## 1. Create a Spectra object with the 3 mzML files from experiment
##    PXD022816.
## 2. Visualise the 3 chromatograms corresponding to the 3
##    acquisitions in that experiment.


px2 <- PXDataset("PXD022816")
mzmls <- head(grep("mzML", pxfiles(px2), value = TRUE), n = 3)
(mzmls <- pxget(px2, mzmls))

(PXD022816 <- Spectra(mzmls))

library(ggplot2)
library(dplyr)

filterMsLevel(PXD022816, 1L) %>%
    spectraData(c("rtime", "totIonCurrent", "dataOrigin")) %>%
    as_tibble() %>%
    mutate(file =  sub("^.+_", "", basename(dataOrigin))) %>%
    ggplot(aes(x = rtime,
               y = totIonCurrent,
               colour = file)) +
    geom_line()
