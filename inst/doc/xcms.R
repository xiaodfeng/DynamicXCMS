## ----biocstyle, echo = FALSE, results = "asis"--------------------------------
BiocStyle::markdown()

## ----init, message = FALSE, echo = FALSE, results = "hide"--------------------
## Silently loading all packages
library(BiocStyle)
library(xcms)
library(faahKO)
library(pander)
## Use socket based parallel processing on Windows systems
## if (.Platform$OS.type == "unix") {
##     register(bpstart(MulticoreParam(3)))
## } else {
##     register(bpstart(SnowParam(3)))
## }
register(SerialParam())


## ----load-libs-pheno, message = FALSE-----------------------------------------
library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)
library(SummarizedExperiment)

## Get the full path to the CDF files
cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)[c(1, 2, 5, 6, 7, 8, 11, 12)]
## Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(cdfs), pattern = ".CDF",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("KO", 4), rep("WT", 4)),
                 stringsAsFactors = FALSE)

## ----load-with-msnbase, message = FALSE---------------------------------------
raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")

## ----subsetting, message = FALSE, echo = TRUE---------------------------------
raw_data <- filterRt(raw_data, c(2500, 3500))

## ----data-inspection-rtime, message = FALSE-----------------------------------
head(rtime(raw_data))

## ----data-inspection-mz, message = FALSE--------------------------------------
mzs <- mz(raw_data)

## Split the list by file
mzs_by_file <- split(mzs, f = fromFile(raw_data))

length(mzs_by_file)

## ----data-inspection-bpc, message = FALSE, fig.align = "center", fig.width = 12, fig.height = 6----
## Get the base peak chromatograms. This reads data from the files.
bpis <- chromatogram(raw_data, aggregationFun = "max")
## Define colors for the two groups
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("KO", "WT")

## Plot all chromatograms.
plot(bpis, col = group_colors[raw_data$sample_group])

## ----data-inspection-chromatogram, message = FALSE----------------------------
bpi_1 <- bpis[1, 1]
head(rtime(bpi_1))
head(intensity(bpi_1))

## ----data-inspection-tic-boxplot, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 4, fig.cap = "Distribution of total ion currents per file."----
## Get the total ion current by file
tc <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tc, col = group_colors[raw_data$sample_group],
        ylab = "intensity", main = "Total ion current")

## ----data-inspection-bpc-heatmap, message = FALSE, fig.align = "center", fig.width = 7, fig.height = 6, fig.cap = "Grouping of samples based on similarity of their base peak chromatogram."----
## Bin the BPC
bpis_bin <- MSnbase::bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = raw_data$sample_group)
rownames(ann) <- raw_data$sample_name

## Perform the cluster analysis
pheatmap(cormat, annotation = ann,
         annotation_color = list(group = group_colors))

## ----peak-detection-plot-eic, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 5, fig.cap = "Extracted ion chromatogram for one peak."----
## Define the rt and m/z range of the peak area
rtr <- c(2700, 2900)
mzr <- c(334.9, 335.1)
## extract the chromatogram
chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
plot(chr_raw, col = group_colors[chr_raw$sample_group])

## ----peak-detection-plot-ms-data, message = FALSE, warning = FALSE, fig.aligh = "center", fig.width = 14, fig.height = 14, fig.cap = "Visualization of the raw MS data for one peak. For each plot: upper panel: chromatogram plotting the intensity values against the retention time, lower panel m/z against retention time plot. The individual data points are colored according to the intensity."----
raw_data %>%
    filterRt(rt = rtr) %>%
    filterMz(mz = mzr) %>%
    plot(type = "XIC")

## ----peak-detection-eic, message = FALSE--------------------------------------
xchr <- findChromPeaks(chr_raw, param = CentWaveParam(snthresh = 2))

## ----peak-detection-eic-chromPeaks--------------------------------------------
head(chromPeaks(xchr))

## ----peak-detection-chromatogram-chromPeakData--------------------------------
chromPeakData(xchr)

## ----peak-detection-eic-plot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Signal for an example peak. Red and blue colors represent KO and wild type samples, respectively. Peak area of identified chromatographic peaks are highlighted in the sample group color."----
sample_colors <- group_colors[xchr$sample_group]
plot(xchr, col = sample_colors,
     peakBg = sample_colors[chromPeaks(xchr)[, "column"]])


## ----peak-detection-centwave, message = FALSE, results = "hide"---------------
cwp <- CentWaveParam(peakwidth = c(20, 80), noise = 5000,
                     prefilter = c(6, 5000))
xdata <- findChromPeaks(raw_data, param = cwp)

## ----peak-detection-chromPeaks, message = FALSE-------------------------------
head(chromPeaks(xdata))

## ----peak-detection-chromPeakData---------------------------------------------
chromPeakData(xdata)

## ----peak-postprocessing, message = FALSE-------------------------------------
mpp <- MergeNeighboringPeaksParam(expandRt = 4)
xdata_pp <- refineChromPeaks(xdata, mpp)

## ----peak-postprocessing-merged, fig.widht = 10, fig.height = 5, fig.cap = "Result from the peak refinement step. Left: data before processing, right: after refinement. The splitted peak was merged into one."----
mzr_1 <- 305.1 + c(-0.01, 0.01)
chr_1 <- chromatogram(filterFile(xdata, 1), mz = mzr_1)
chr_2 <- chromatogram(filterFile(xdata_pp, 1), mz = mzr_1)
par(mfrow = c(1, 2))
plot(chr_1)
plot(chr_2)

## ----peak-postprocessing-not-merged, fig.widht = 10, fig.height = 5, fig.cap = "Result from the peak refinement step. Left: data before processing, right: after refinement. The peaks were not merged."----
mzr_1 <- 496.2 + c(-0.01, 0.01)
chr_1 <- chromatogram(filterFile(xdata, 1), mz = mzr_1)
chr_2 <- chromatogram(filterFile(xdata_pp, 1), mz = mzr_1)
par(mfrow = c(1, 2))
plot(chr_1)
plot(chr_2)

## ----peak-postprocessing-chr, fig..width = 5, fig.height = 5------------------
res <- refineChromPeaks(chr_1, MergeNeighboringPeaksParam(minProp = 0.05))
chromPeaks(res)
plot(res)

## -----------------------------------------------------------------------------
xdata <- xdata_pp

## ----peak-detection-peaks-per-sample, message = FALSE, results = "asis"-------
summary_fun <- function(z)
    c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))

T <- lapply(split.data.frame(
    chromPeaks(xdata), f = chromPeaks(xdata)[, "sample"]),
    FUN = summary_fun)
T <- do.call(rbind, T)
rownames(T) <- basename(fileNames(xdata))
pandoc.table(
    T,
    caption = paste0("Summary statistics on identified chromatographic",
                     " peaks. Shown are number of identified peaks per",
                     " sample and widths/duration of chromatographic ",
                     "peaks."))

## ----peak-detection-chrom-peaks-plot, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 8, fig.cap = "Identified chromatographic peaks in the m/z by retention time space for one sample."----
plotChromPeaks(xdata, file = 3)

## ----peak-detection-chrom-peak-image, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Frequency of identified chromatographic peaks along the retention time axis. The frequency is color coded with higher frequency being represented by yellow-white. Each line shows the peak frequency for one file."----
plotChromPeakImage(xdata)

## ----peak-detection-eic-example-peak, message = FALSE-------------------------
chr_ex <- chromatogram(xdata, mz = mzr, rt = rtr)
chromPeaks(chr_ex)

## ----peak-detection-highlight-chrom-peaks-plot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Signal for an example peak. Red and blue colors represent KO and wild type samples, respectively. The rectangles indicate the identified chromatographic peaks per sample."----
sample_colors <- group_colors[chr_ex$sample_group]
plot(chr_ex, col = sample_colors, peakType = "rectangle",
     peakCol = sample_colors[chromPeaks(chr_ex)[, "sample"]],
     peakBg = NA)

## ----peak-detection-highlight-chrom-peaks-plot-polygon, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Signal for an example peak. Red and blue colors represent KO and wild type samples, respectively. The signal area of identified chromatographic peaks are filled with a color."----
plot(chr_ex, col = group_colors[chr_raw$sample_group], lwd = 2,
     peakBg = sample_colors[chromPeaks(chr_ex)[, "sample"]])

## ----peak-detection-chrom-peak-table-selected, message = FALSE, results = "asis"----
pander(chromPeaks(xdata, mz = mzr, rt = rtr),
       caption = paste("Identified chromatographic peaks in a selected ",
                       "m/z and retention time range."))

## ----peak-detection-chrom-peak-intensity-boxplot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Peak intensity distribution per sample."----
## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xdata)[, "into"]),
              f = chromPeaks(xdata)[, "sample"])
boxplot(ints, varwidth = TRUE, col = group_colors[xdata$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities")
grid(nx = NA, ny = NULL)

## ----alignment-obiwarp, message = FALSE, results = "hide"---------------------
xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

## ----alignment-rtime, message = FALSE-----------------------------------------
## Extract adjusted retention times
head(adjustedRtime(xdata))

## Or simply use the rtime method
head(rtime(xdata))

## ----alignment-obiwarp-plot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Obiwarp aligned data. Base peak chromatogram after alignment (top) and difference between adjusted and raw retention times along the retention time axis (bottom)."----
## Get the base peak chromatograms.
bpis_adj <- chromatogram(xdata, aggregationFun = "max", include = "none")
par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(bpis_adj, col = group_colors[bpis_adj$sample_group])
## Plot also the difference of adjusted to raw retention time.
plotAdjustedRtime(xdata, col = group_colors[xdata$sample_group])

## ----alignment-peak-groups-example-peak, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Example extracted ion chromatogram before (top) and after alignment (bottom)."----
par(mfrow = c(2, 1))
## Plot the raw data
plot(chr_raw, col = group_colors[chr_raw$sample_group])

## Extract the chromatogram from the adjusted object
chr_adj <- chromatogram(xdata, rt = rtr, mz = mzr)
plot(chr_adj, col = group_colors[chr_raw$sample_group], peakType = "none")

## ----subset-define, message = FALSE, warning = FALSE--------------------------
xdata <- dropAdjustedRtime(xdata)

## Define the experimental layout
xdata$sample_type <- "study"
xdata$sample_type[c(1, 4, 7)] <- "QC"

## ----alignment-subset, message = FALSE, warning = FALSE-----------------------
## Initial peak grouping. Use sample_type as grouping variable
pdp_subs <- PeakDensityParam(sampleGroups = xdata$sample_type,
                             minFraction = 0.9)
xdata <- groupChromPeaks(xdata, param = pdp_subs)

## Define subset-alignment options and perform the alignment
pgp_subs <- PeakGroupsParam(minFraction = 0.85,
                            subset = which(xdata$sample_type == "QC"),
                            subsetAdjust = "average", span = 0.4)
xdata <- adjustRtime(xdata, param = pgp_subs)

## ----alignment-subset-plot-2, message = FALSE, warning = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Subset-alignment results with option average. Difference between adjusted and raw retention times along the retention time axis. Samples on which the alignment models were estimated are shown in green, study samples in grey."----
clrs <- rep("#00000040", 8)
clrs[xdata$sample_type == "QC"] <- c("#00ce0080")
par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 0.5))
plot(chromatogram(xdata, aggregationFun = "sum"),
     col = clrs, peakType = "none")
plotAdjustedRtime(xdata, col = clrs, peakGroupsPch = 1,
                  peakGroupsCol = "#00ce0040")

## ----correspondence-example-slice, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Example for peak density correspondence. Upper panel: chromatogram for an mz slice with multiple chromatographic peaks. lower panel: identified chromatographic peaks at their retention time (x-axis) and index within samples of the experiments (y-axis) for different values of the bw parameter. The black line represents the peak density estimate. Grouping of peaks (based on the provided settings) is indicated by grey rectangles."----
## Define the mz slice.
mzr <- c(305.05, 305.15)

## Extract and plot the chromatograms
chr_mzr <- chromatogram(xdata, mz = mzr)
## Define the parameters for the peak density method
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 30)
plotChromPeakDensity(chr_mzr, col = sample_colors, param = pdp,
                     peakBg = sample_colors[chromPeaks(chr_mzr)[, "sample"]],
                     peakCol = sample_colors[chromPeaks(chr_mzr)[, "sample"]],
                     peakPch = 16)


## ----correspondence, message = FALSE------------------------------------------
## Perform the correspondence
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.4, bw = 30)
xdata <- groupChromPeaks(xdata, param = pdp)

## ----quantify-----------------------------------------------------------------
res <- quantify(xdata, value = "into")

## ----quantify-colData---------------------------------------------------------
colData(res)

## ----quantify-rowData---------------------------------------------------------
rowData(res)

## -----------------------------------------------------------------------------
head(assay(res))

## ----correspondence-featureDefs, message = FALSE------------------------------
## Extract the feature definitions
featureDefinitions(xdata)

## ----correspondence-feature-values, message = FALSE---------------------------
## Extract the into column for each feature.
head(featureValues(xdata, value = "into"))

## ----fill-chrom-peaks, message = FALSE----------------------------------------
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam())

head(featureValues(xdata))

## ----fill-chrom-peaks-compare, message = FALSE--------------------------------
## Missing values before filling in peaks
apply(featureValues(xdata, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))

## Missing values after filling in peaks
apply(featureValues(xdata), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))


## ----export-result, eval = FALSE, echo = FALSE--------------------------------
#  save(xdata, file = "xdata.RData")

## ----featureSummary, message = FALSE------------------------------------------
head(featureSummary(xdata, group = xdata$sample_group))

## -----------------------------------------------------------------------------
assays(res)$raw_filled <- featureValues(xdata, filled = TRUE)

## -----------------------------------------------------------------------------
assayNames(res)

head(assay(res, "raw"))
head(assay(res, "raw_filled"))

## ----featureChromatograms, message = FALSE------------------------------------
feature_chroms <- featureChromatograms(xdata, features = 1:4)

feature_chroms

## ----feature-eic, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 8, fig.cap = "Extracted ion chromatograms for features 1 to 4."----
plot(feature_chroms, col = sample_colors,
     peakBg = sample_colors[chromPeaks(feature_chroms)[, "sample"]])


## -----------------------------------------------------------------------------
eic_2 <- feature_chroms[2, ]
chromPeaks(eic_2)

## ----final-pca, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 8, fig.cap = "PCA for the faahKO data set, un-normalized intensities."----
## Extract the features and log2 transform them
ft_ints <- log2(assay(res, "raw_filled"))

## Perform the PCA omitting all features with an NA in any of the
## samples. Also, the intensities are mean centered.
pc <- prcomp(t(na.omit(ft_ints)), center = TRUE)

## Plot the PCA
cols <- group_colors[xdata$sample_group]
pcSummary <- summary(pc)
plot(pc$x[, 1], pc$x[,2], pch = 21, main = "",
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = "darkgrey", bg = cols, cex = 2)
grid()
text(pc$x[, 1], pc$x[,2], labels = xdata$sample_name, col = "darkgrey",
     pos = 3, cex = 2)


## ----processhistory, message = FALSE------------------------------------------
processHistory(xdata)

## ----processhistory-select, message = FALSE-----------------------------------
ph <- processHistory(xdata, type = "Retention time correction")

ph

## ----processhistory-param, message = FALSE------------------------------------
## Access the parameter
processParam(ph[[1]])


## ----subset-filterFile, message = FALSE---------------------------------------
subs <- filterFile(xdata, file = c(2, 4))

## Do we have identified chromatographic peaks?
hasChromPeaks(subs)

## ----subset-filterFile-2, message = FALSE-------------------------------------
## Do we still have features?
hasFeatures(subs)

## Do we still have adjusted retention times?
hasAdjustedRtime(subs)

## ----subset-filterFile-3, message = FALSE-------------------------------------
subs <- filterFile(xdata, keepAdjustedRtime = TRUE)

hasAdjustedRtime(subs)

## ----subset-filterRt, message = FALSE-----------------------------------------
subs <- filterRt(xdata, rt = c(3000, 3500))

range(rtime(subs))

## ----subset-filterRt-2, message = FALSE---------------------------------------
hasAdjustedRtime(subs)

## ----subset-filterRt-3, message = FALSE---------------------------------------
hasChromPeaks(subs)

range(chromPeaks(subs)[, "rt"])

## ----subset-bracket, message = FALSE, warning = FALSE-------------------------
## Extract all data from the 3rd file.
one_file <- filterFile(xdata, file = 3)

one_file_2 <- xdata[fromFile(xdata) == 3]

## Is the content the same?
all.equal(one_file[[2]], one_file_2[[2]])

## ----subset-bracket-peaks, message = FALSE------------------------------------
## Subsetting with filterFile preserves chromatographic peaks
head(chromPeaks(one_file))

## Subsetting with [ not
head(chromPeaks(one_file_2))

## ----subset-bracket-keepAdjustedRtime, message = FALSE, warnings = FALSE------
subs <- xdata[20:30, keepAdjustedRtime = TRUE]

hasAdjustedRtime(subs)

## Access adjusted retention times:
rtime(subs)

## Access raw retention times:
rtime(subs, adjusted = FALSE)

## ----subset-double-bracket, message = FALSE-----------------------------------
## Extract a single spectrum
xdata[[14]]

## ----subset-split, message = FALSE--------------------------------------------
x_list <- split(xdata, f = fromFile(xdata), keepAdjustedRtime = TRUE)

lengths(x_list)

lapply(x_list, hasAdjustedRtime)

## ----multicore, message = FALSE, eval = FALSE---------------------------------
#  register(bpstart(MulticoreParam(2)))

## ----snow, message = FALSE, eval = FALSE--------------------------------------
#  register(bpstart(SnowParam(2)))

