.libPaths("d:/PackageR/4.1/") # change the default package directory
library(devtools)
install("d:/github/DynamicXCMS",update=F)
# install("d:/PackageR/xcms/xcms_3.16.1/",dep=TRUE)

library(xcms)

