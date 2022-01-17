## ----biocstyle, echo = FALSE, results = "asis"--------------------------------
BiocStyle::markdown()

## ----init, message = FALSE, echo = FALSE, results = "hide"--------------------
## Silently loading all packages
library(BiocStyle)
library(xcms)
library(Spectra)
library(pander)
register(SerialParam())


## ----load-dda-data, message = FALSE-------------------------------------------
library(xcms)

dda_file <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML",
                        package = "msdata")
dda_data <- readMSData(dda_file, mode = "onDisk")

## ----subset-dda, echo = FALSE, message = FALSE, eval = TRUE-------------------
#' Silently sub-setting the object to speed-up analysis
dda_data <- filterRt(dda_data, rt = c(200, 600))

## ----dda-table-mslevel--------------------------------------------------------
table(msLevel(dda_data))

## ----precursor----------------------------------------------------------------
library(magrittr)

dda_data %>%
    filterMsLevel(2L) %>%
    precursorMz() %>%
    head()

## ----precursor-intensity------------------------------------------------------
dda_data %>%
    filterMsLevel(2L) %>%
    precursorIntensity() %>%
    head()

## ----estimate-precursor-------------------------------------------------------
prec_int <- xcms::estimatePrecursorIntensity(dda_data)

## ----set-precursor-intensity--------------------------------------------------
fData(dda_data)$precursorIntensity <- prec_int

dda_data %>%
    filterMsLevel(2L) %>%
    precursorIntensity() %>%
    head()

## ----dda-find-chrom-peaks-ms1, message = FALSE--------------------------------
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                     peakwidth = c(3, 30))
dda_data <- findChromPeaks(dda_data, param = cwp)

## ----dda-spectra, message = FALSE, eval = TRUE--------------------------------
library(Spectra)
dda_spectra <- chromPeakSpectra(
    dda_data, msLevel = 2L, return.type = "Spectra")
dda_spectra

## ----peak_id, eval = TRUE-----------------------------------------------------
dda_spectra$peak_id

## ----dda-ms2-example, message = FALSE, eval = TRUE----------------------------
ex_mz <- 304.1131
chromPeaks(dda_data, mz = ex_mz, ppm = 20)

## ----dda-ms2-get-ms2, message = FALSE, eval = TRUE----------------------------
ex_id <- rownames(chromPeaks(dda_data, mz = ex_mz, ppm = 20))
ex_spectra <- dda_spectra[dda_spectra$peak_id == ex_id]
ex_spectra

## ----dda-ms2-consensus, message = FALSE, eval = TRUE--------------------------
ex_spectrum <- combineSpectra(ex_spectra, FUN = combinePeaks, ppm = 20,
                              peaks = "intersect", minProp = 0.8,
                              intensityFun = median, mzFun = median,
                              f = ex_spectra$peak_id)
ex_spectrum

## ----dda-ms2-consensus-plot, message = FALSE, fig.cap = "Consensus MS2 spectrum created from all measured MS2 spectra for ions of chromatographic peak CP53.", fig.width = 8, fig.height = 8, eval = TRUE----
plotSpectra(ex_spectrum)

## ----normalize, eval = TRUE---------------------------------------------------
norm_fun <- function(z, ...) {
    z[, "intensity"] <- z[, "intensity"] /
        max(z[, "intensity"], na.rm = TRUE) * 100
    z
}
ex_spectrum <- addProcessing(ex_spectrum, FUN = norm_fun)

## ----dda-ms2-metlin-match, fig.cap = "Mirror plots for the candidate MS2 spectrum against Flumanezil (left) and Fenamiphos (right). The upper panel represents the candidate MS2 spectrum, the lower the target MS2 spectrum. Matching peaks are indicated with a dot.", fig.width = 12, fig.height = 6, eval = TRUE----
library(MsBackendMgf)
flumanezil <- Spectra(
    system.file("mgf", "metlin-2724.mgf", package = "xcms"),
    source = MsBackendMgf())
fenamiphos <- Spectra(
    system.file("mgf", "metlin-72445.mgf", package = "xcms"),
    source = MsBackendMgf())

par(mfrow = c(1, 2))
plotSpectraMirror(ex_spectrum, flumanezil[3], main = "against Flumanezil",
                  ppm = 40)
plotSpectraMirror(ex_spectrum, fenamiphos[3], main = "against Fenamiphos",
                  ppm = 40)

## ----dda-ms2-dotproduct, eval = TRUE------------------------------------------
compareSpectra(ex_spectrum, flumanezil, ppm = 40)
compareSpectra(ex_spectrum, fenamiphos, ppm = 40)

## ----load-swath-data, message = FALSE-----------------------------------------
swath_file <- system.file("TripleTOF-SWATH",
                          "PestMix1_SWATH.mzML",
                          package = "msdata")

swath_data <- readMSData(swath_file, mode = "onDisk")

## ----echo = FALSE, message = FALSE--------------------------------------------
swath_data <- filterRt(swath_data, rt = c(200, 600))

## ----swath-table-mslevel------------------------------------------------------
table(msLevel(swath_data))

## ----fdata-isolationwindow----------------------------------------------------
head(fData(swath_data)[, c("isolationWindowTargetMZ",
                           "isolationWindowLowerOffset",
                           "isolationWindowUpperOffset",
                           "msLevel", "retentionTime")])

head(isolationWindowLowerMz(swath_data))
head(isolationWindowUpperMz(swath_data))

## -----------------------------------------------------------------------------
table(isolationWindowTargetMz(swath_data))

## ----find-chrom-peaks-ms1, message = FALSE------------------------------------
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaks(swath_data, param = cwp)

## ----find-chrom-peaks-ms2, message = FALSE------------------------------------
cwp <- CentWaveParam(snthresh = 3, noise = 10, ppm = 10,
                     peakwidth = c(3, 30))
swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp)

## -----------------------------------------------------------------------------
chromPeakData(swath_data)

## -----------------------------------------------------------------------------
table(chromPeakData(swath_data)$isolationWindow)

## ----fena-extract-peak--------------------------------------------------------
fenamiphos_mz <- 304.113077
fenamiphos_ms1_peak <- chromPeaks(swath_data, mz = fenamiphos_mz, ppm = 2)
fenamiphos_ms1_peak

## ----fena-identify-ms2--------------------------------------------------------
keep <- chromPeakData(swath_data)$isolationWindowLowerMz < fenamiphos_mz &
        chromPeakData(swath_data)$isolationWindowUpperMz > fenamiphos_mz

## ----fena-check-rt------------------------------------------------------------
keep <- keep &
    chromPeaks(swath_data)[, "rtmin"] < fenamiphos_ms1_peak[, "rt"] &
    chromPeaks(swath_data)[, "rtmax"] > fenamiphos_ms1_peak[, "rt"]

fenamiphos_ms2_peak <- chromPeaks(swath_data)[which(keep), ]

## ----fena-eic-extract, warning = FALSE----------------------------------------
rtr <- fenamiphos_ms1_peak[, c("rtmin", "rtmax")]
mzr <- fenamiphos_ms1_peak[, c("mzmin", "mzmax")]
fenamiphos_ms1_chr <- chromatogram(swath_data, rt = rtr, mz = mzr)

rtr <- fenamiphos_ms2_peak[, c("rtmin", "rtmax")]
mzr <- fenamiphos_ms2_peak[, c("mzmin", "mzmax")]
fenamiphos_ms2_chr <- chromatogram(
    filterIsolationWindow(swath_data, mz = fenamiphos_mz),
    rt = rtr, mz = mzr, msLevel = 2L)

## ----fena-eic-plot, fig.width = 10, fig.height = 5, fig.cap = "Extracted ion chromatograms for Fenamiphos from MS1 (blue) and potentially related signal in MS2 (grey)."----
plot(rtime(fenamiphos_ms1_chr[1, 1]),
     intensity(fenamiphos_ms1_chr[1, 1]),
     xlab = "retention time [s]", ylab = "intensity", pch = 16,
     ylim = c(0, 5000), col = "blue", type = "b", lwd = 2)
#' Add data from all MS2 peaks
tmp <- lapply(fenamiphos_ms2_chr@.Data,
              function(z) points(rtime(z), intensity(z),
                                 col = "#00000080",
                                 type = "b", pch = 16))

## ----fena-cor-----------------------------------------------------------------
correlate(fenamiphos_ms2_chr[1, 1],
          fenamiphos_ms1_chr[1, 1], align = "approx")

## ----reconstruct-ms2, message = FALSE-----------------------------------------
swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.9,
                                             return.type = "Spectra")
swath_spectra

## -----------------------------------------------------------------------------
swath_spectra$ms2_peak_id
swath_spectra$peak_id

## ----fena-swath-peak----------------------------------------------------------
fenamiphos_swath_spectrum <- swath_spectra[
    swath_spectra$peak_id == rownames(fenamiphos_ms1_peak)]

## -----------------------------------------------------------------------------
fenamiphos_swath_spectrum <- addProcessing(fenamiphos_swath_spectrum,
                                           norm_fun)

## ----fena-swath-plot, fig.cap = "Mirror plot comparing the reconstructed MS2 spectrum for Fenamiphos (upper panel) against the measured spectrum from the DDA data and the Fenamiphhos spectrum from Metlin.", fig.width = 12, fig.height = 6----
par(mfrow = c(1, 2))
plotSpectraMirror(fenamiphos_swath_spectrum, ex_spectrum,
     ppm = 50, main = "against DDA")
plotSpectraMirror(fenamiphos_swath_spectrum, fenamiphos[2],
     ppm = 50, main = "against Metlin")

## -----------------------------------------------------------------------------
pk_ids <- fenamiphos_swath_spectrum$ms2_peak_id[[1]]
pk_ids

## -----------------------------------------------------------------------------
rt_range <- chromPeaks(swath_data)[pk_ids, c("rtmin", "rtmax")]
mz_range <- chromPeaks(swath_data)[pk_ids, c("mzmin", "mzmax")]

pmz <- precursorMz(fenamiphos_swath_spectrum)[1]
swath_data_iwindow <- filterIsolationWindow(swath_data, mz = pmz)
ms2_eics <- chromatogram(swath_data_iwindow, rt = rt_range,
                         mz = mz_range, msLevel = 2L)

## ----pro-swath----------------------------------------------------------------
prochloraz_mz <- 376.0381

prochloraz_ms1_peak <- chromPeaks(swath_data, msLevel = 1L,
                                  mz = prochloraz_mz, ppm = 5)
prochloraz_ms1_peak

prochloraz_swath_spectrum <- swath_spectra[
    swath_spectra$peak_id == rownames(prochloraz_ms1_peak)]

## ----pro-dda------------------------------------------------------------------
prochloraz_dda_peak <- chromPeaks(dda_data, msLevel = 1L,
                                  mz = prochloraz_mz, ppm = 5)
prochloraz_dda_peak

## ----pro-dda-ms2--------------------------------------------------------------
prochloraz_dda_spectra <- dda_spectra[
    dda_spectra$peak_id == rownames(prochloraz_dda_peak)]
prochloraz_dda_spectra

## ----pro-dda-consensus--------------------------------------------------------
prochloraz_dda_spectrum <- combineSpectra(
    prochloraz_dda_spectra, FUN = combinePeaks, ppm = 20,
    peaks = "intersect", minProp = 0.8, intensityFun = median, mzFun = median,
    f = prochloraz_dda_spectra$peak_id)

## ----prochloraz-metlin--------------------------------------------------------
prochloraz <- Spectra(
    system.file("mgf", "metlin-68898.mgf", package = "xcms"),
    source = MsBackendMgf())

## ----pro-swath-plot, fig.cap = "Mirror plot comparing the reconstructed MS2 spectrum for Prochloraz (upper panel) against the measured spectrum from the DDA data and the Prochloraz spectrum from Metlin.", fig.width = 12, fig.height = 6----
prochloraz_swath_spectrum <- addProcessing(prochloraz_swath_spectrum, norm_fun)
prochloraz_dda_spectrum <- addProcessing(prochloraz_dda_spectrum, norm_fun)

par(mfrow = c(1, 2))
plotSpectraMirror(prochloraz_swath_spectrum, prochloraz_dda_spectrum,
                  ppm = 40, main = "against DDA")
plotSpectraMirror(prochloraz_swath_spectrum, prochloraz[2],
                  ppm = 40, main = "against Metlin")

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

