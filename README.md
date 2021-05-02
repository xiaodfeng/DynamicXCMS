[![Build Status](https://travis-ci.com/xiaodfeng/DynamicXCMS.svg?branch=main)](https://travis-ci.com/xiaodfeng/DynamicXCMS)

# The `dynamic xcms` package (based on xcms 3.8.2)

<img align = "right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/xcms/xcms.png" height="200">


## How to track the code changes.
The original source codes of xcms 3.8.2 have been downloaded and uploaded to https://github.com/xiaodfeng/DynamicXCMS, with the commit “3.8.2, original xcms codes”. Based on this stable version, we implement the dynamic theory into the package, with the commit “3.8.2.1, dynamic xcms”. In this dynamic xcms package, seven files were changed: DESCRIPTION, DynamicXCMS.Rproj, DataClasses.R, do_findChromPeaks-functions.R, functions-Params.R, methods-Params.R, and mzROI.c, detailed changes can be accessed in the depository. 
## How to install the package.
Download the source codes from https://github.com/xiaodfeng/DynamicXCMS or from Appendix C of the paper.
Install the downloaded codes with the following commands:
*library(devtools)*
*install("../DynamicXCMS ") # install the package locally, not from the website. See the output below*
*library(xcms)*

## How to use the package.
In the original xcms, the peak detection mass tolerance (PDMT) was set by the parameter ppm. While in the dynamic xcms, the PDMT was set by parameters A, ppm, and instrument. The A in dynamic xcms can be calculated by the reference m/z and the reference mass resolving power of the .raw data and equals 4.2897〖∙10〗^(-7) according to equation 20 of the manuscript. The ppm in dynamic xcms was used to set the mass fluctuation (MF) and equals 1 according to equation 20 of the manuscript. The instrument in dynamic xcms indicates the instrument types, with FTICR=1, Orbitrap=2, Q-TOF=3, and Quadrupole=4. With this research's Orbitrap dataset, the parameters for peak detection were set below:
*c*wp <- CentWaveParam(*
  *A = 4.289723e-07, ppm=1,  Instrument=2,*
  *peakwidth=c(2.4,30), snthresh = 10, noise=100000,* 
  prefilter=c(1, 200000),  firstBaselineCheck = TRUE, integrate=2)*
Peak detection with dynamic xcms, see output below:
*xdatad <- findChromPeaks(Raw_data, param = cwp,BPPARAM = BPPARAM)*


## Discussions and suggestions are welcome:
https://github.com/xiaodfeng/DynamicXCMS/issues

