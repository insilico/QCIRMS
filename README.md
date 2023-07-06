# QCIRMS
Quality Control Analysis for Isotope Ratio Mass Spectrometry Data

### To install:

>library(devtools)
>install_github("insilico/QCIRMS")  
>library(QCIRMS)

### Dependencies
```
install.packages(c('dplyr', 'pracma'))
```

```
library(devtools)
install_github('isoverse/isoreader')
```

### Examples

```
# example dxf data file provided in package: inst/extdata
library(QCIRMS)
library(isoreader)
dxf_file <- "170506_NaHCO3 L + NaCl U_.dxf"
dataPath <- system.file("extdata/dxf_files/abiotic", dxf_file, package = "QCIRMS")
file.summ<-read_summary(dataPath)
```
### Abstract

#### Contact
[lily-clough@utulsa.edu](lily-clough@utulsa.edu)
