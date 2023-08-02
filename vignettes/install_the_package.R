#### Xiaodong 202308 ####
# This script is used for installing the package, as the latest version of 
# MSnbase 2.27.1 (https://github.com/lgatto/MSnbase/blob/master/DESCRIPTION) 
# has changed the the Spectra and Chromatograms class into MSpectra and MChromatograms class
# Thus installing the XCMS 3.8.2.1 based on MSnbase 2.27.1 will result in an error. 
# We need to stick to the MSnbase-2.14.2 which can be downloaded via https://code.bioconductor.org/browse/MSnbase/blob/RELEASE_3_11/DESCRIPTION
# In the long run, the DynamicXCMS needs to be implanted based on the latest version of XCMS version 4
library(devtools)
devtools::install("d:/OneDrive - NKI/github/PackageR/MSnbase-2.14.2",update = FALSE)
devtools::install("d:/OneDrive - NKI/github/DynamicXCMS",update = FALSE)
library(MSnbase)
library(xcms)
sessionInfo()
