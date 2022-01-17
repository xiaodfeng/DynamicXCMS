## ----biocstyle, echo = FALSE, results = "asis"--------------------------------
BiocStyle::markdown()
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## ----init, results = "hide", echo = FALSE-------------------------------------
## Silently loading all packages
library(BiocStyle)
library(xcms)
library(MsFeatures)
register(SerialParam())


## ----load-data----------------------------------------------------------------
library(xcms)
library(faahKO)
library(MsFeatures)

data("xdata")
## Update the path to the files for the local system
dirname(xdata) <- c(rep(system.file("cdf", "KO", package = "faahKO"), 4),
                    rep(system.file("cdf", "WT", package = "faahKO"), 4))

## ----fdev---------------------------------------------------------------------
featureDefinitions(xdata)

## ----filled-not-filled--------------------------------------------------------
head(featureValues(xdata, filled = FALSE))
head(featureValues(xdata, filled = TRUE))

## ----feature-rt-mz-plot, fig.width = 8, fig.height = 6, fig.cap = "Plot of retention times and m/z for all features in the data set."----
plot(featureDefinitions(xdata)$rtmed, featureDefinitions(xdata)$mzmed,
     xlab = "retention time", ylab = "m/z", main = "features",
     col = "#00000080")
grid()

## -----------------------------------------------------------------------------
xdata <- groupFeatures(xdata, param = SimilarRtimeParam(10))

## -----------------------------------------------------------------------------
table(featureGroups(xdata))

## ----feature-groups-rtime-plot, fig.width = 8, fig.height = 6, fig.cap = "Feature groups defined with a rt window of 10 seconds"----
plotFeatureGroups(xdata)
grid()

## ----repeat-------------------------------------------------------------------
## Remove previous feature grouping results to repeat the rtime-based
## feature grouping with different setting
featureDefinitions(xdata)$feature_group <- NULL

## Repeat the grouping
xdata <- groupFeatures(xdata, SimilarRtimeParam(20))
table(featureGroups(xdata))

## ----feature-groups-rtime-plot2, fig.width = 8, fig.height = 6, fig.cap = "Feature groups defined with a rt window of 20 seconds"----
plotFeatureGroups(xdata)
grid()

## ----abundance-correlation-heatmap, fig.cap = "Correlation of features based on feature abundances.", fig.width = 6, fig.height = 16----
library(pheatmap)
fvals <- log2(featureValues(xdata, filled = TRUE))

cormat <- cor(t(fvals), use = "pairwise.complete.obs")
ann <- data.frame(fgroup = featureGroups(xdata))
rownames(ann) <- rownames(cormat)

res <- pheatmap(cormat, annotation_row = ann, cluster_rows = TRUE,
                cluster_cols = TRUE)

## ----abundance-correlation----------------------------------------------------
xdata <- groupFeatures(xdata, AbundanceSimilarityParam(threshold = 0.7,
                                                       transform = log2),
                       filled = TRUE)
table(featureGroups(xdata))

## ----abundance-correlation-fg040, fig.width = 8, fig.height = 8, fig.cap = "Pairwise correlation plot for all features initially grouped into the feature group FG.040."----
fts <- grep("FG.040", featureGroups(xdata))
pairs(t(fvals[fts, ]), gap = 0.1, main = "FG.040")

## ----correlate-eic, message = FALSE-------------------------------------------
xdata <- groupFeatures(xdata, EicSimilarityParam(threshold = 0.7, n = 2))

## ----correlate-eic-result-----------------------------------------------------
table(featureGroups(xdata))

## -----------------------------------------------------------------------------
fts <- grep("FG.008.001", featureGroups(xdata))
eics <- featureChromatograms(xdata, features = fts,
                             filled = TRUE, n = 1)

## ----example-1-eic, fig.width = 8, fig.height = 6, fig.cap = "EICs of features from feature group FG.008.001 in the same sample. Shown are the actual intensities (left) and intensities normalized to 1 (right). Features being part of the same feature group after grouping by EIC similarity are shown in the same color."----
cols <- c("#ff000080", "#00ff0080")
names(cols) <- unique(featureGroups(xdata)[fts])

par(mfrow = c(1, 2))
plotChromatogramsOverlay(eics, col = cols[featureGroups(xdata)[fts]],
                         lwd = 2, peakType = "none")
plotChromatogramsOverlay(normalize(eics),
                         col = cols[featureGroups(xdata)[fts]],
                         lwd = 2, peakType = "none")

## -----------------------------------------------------------------------------
fts <- grep("FG.068.001", featureGroups(xdata))
eics <- featureChromatograms(xdata, features = fts,
                             filled = TRUE, n = 1)

## ----example-2-eic, fig.width = 8, fig.height = 6, fig.cap = "EICs for features from feature group FG.068.001 in the same sample. Shown are the actual intensities (left) and intensities normalized to 1 (right). Features being part of the same feature group after grouping by EIC similarity are shown in the same color."----
cols <- c("#ff000080", "#00ff0080")
names(cols) <- unique(featureGroups(xdata)[fts])

par(mfrow = c(1, 2))
plotChromatogramsOverlay(eics, col = cols[featureGroups(xdata)[fts]],
                         lwd = 2, peakType = "none")
plotChromatogramsOverlay(normalize(eics),
                         col = cols[featureGroups(xdata)[fts]],
                         lwd = 2, peakType = "none")


## ----reset-feature-groups-----------------------------------------------------
featureDefinitions(xdata)$feature_group <- NA_character_

set.seed(123)
fts_idx <- sample(1:nrow(featureDefinitions(xdata)), 30)
featureDefinitions(xdata)$feature_group[fts_idx] <- "FG"

## ----rtime-grouping-----------------------------------------------------------
xdata <- groupFeatures(xdata, SimilarRtimeParam(diffRt = 20))
xdata <- groupFeatures(xdata, AbundanceSimilarityParam(threshold = 0.7))
table(featureGroups(xdata))

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

