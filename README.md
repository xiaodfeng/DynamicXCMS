[![Build Status](https://travis-ci.com/xiaodfeng/DynamicXCMS.svg?branch=main)](https://travis-ci.com/xiaodfeng/DynamicXCMS)

# The dynamic xcms package (based on xcms 3.8.2)

<img align = "right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/xcms/xcms.png" height="200">

## How to track the code changes

The dynamic xcms package (verson 3.8.2.1) is developted based on the stable version of [xcms version 3.8.2 ](https://code.bioconductor.org/browse/xcms/blob/RELEASE_3_10/DESCRIPTION). with the dynamic theory implanted. In this dynamic xcms package, seven files were changed: DESCRIPTION, DynamicXCMS.Rproj, DataClasses.R, do_findChromPeaks-functions.R, functions-Params.R, methods-Params.R, and mzROI.c, detailed changes can be accessed in the commits of the depository. 

## How to install the package.

Download the source codes from https://github.com/xiaodfeng/DynamicXCMS or from Appendix C of PMID: *34172146* into your directory

As the latest version of [MSnbase 2.27.1](https://github.com/lgatto/MSnbase/blob/master/DESCRIPTION) has changed the the *Spectra* and *Chromatograms* class into *MSpectra* and *MChromatograms* class. Thus installing the XCMS 3.8.2.1 based on MSnbase 2.27.1 will result in an error. We need to stick to the [MSnbase-2.14.2](https://code.bioconductor.org/browse/MSnbase/blob/RELEASE_3_11/DESCRIPTION). 

In the long run, the DynamicXCMS needs to be implanted based on the latest version of XCMS version 4.

> library(devtools)
> 
> 
> devtools::install("d:/your directory/MSnbase-2.14.2",update = FALSE)
> 
> 
> devtools::install("d:/your directory/DynamicXCMS",update = FALSE)
> 
> 
> library(MSnbase)
> 
> 
> library(xcms)

## How to use the package.

In the original xcms, the peak detection mass tolerance (PDMT) was set by the parameter ppm. While in the dynamic xcms, the PDMT was set by parameters A, ppm, and instrument. 

- The A in dynamic xcms can be calculated by the reference m/z and the reference mass resolving power of the .raw data. The default value of A is set to $4.2897/10^7$ according to equation (20) of the manuscript PMID: *34172146*. 

- The ppm in dynamic xcms was used to set the mass fluctuation (MF) and the default vaue is set to 1 according to equation (20) of the manuscript PMID: *34172146*.

- The instrument in dynamic xcms indicates the instrument types, with FTICR=1, Orbitrap=2, Q-TOF=3, and Quadrupole=4. 

- With the Orbitrap dataset used in the manuscript, the parameters for peak detection were set below:
  
  > #set the parameters
  > 
  > cwp <- CentWaveParam(A = 4.289723e-07, ppm=1, Instrument=2,peakwidth=c(2.4,30), snthresh = 10, noise=100000, prefilter=c(1, 200000),  firstBaselineCheck = TRUE, integrate=2) 
  > 
  > #Peak detection with dynamic xcms
  > 
  > xdatad <- findChromPeaks(Raw_data, param = cwp,BPPARAM = BPPARAM) 

## Cite

Feng, X., Zhang, W., Kuipers, F., Kema, I., Barcaru, A., & Horvatovich, P. (2021). Dynamic binning peak detection and assessment of various lipidomics liquid chromatography-mass spectrometry pre-processing platforms. *Analytica Chimica Acta*, *31*(0), 1–17. https://doi.org/10.1016/j.aca.2021.338674

## Discussions and suggestions are welcome

https://github.com/xiaodfeng/DynamicXCMS/issues
