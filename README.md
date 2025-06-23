# QCIRMS
Quality Control Analysis for Isotope Ratio Mass Spectrometry Data

### To install:

```
install.packages(c('dplyr', 'pracma'))
install.packages("devtools")
library(devtools)
install_github("isoverse/isoreader")
install_github("insilico/QCIRMS") 
library(QCIRMS)
```

### Dependencies

### Examples

```
# example dxf data file provided in package: inst/extdata
library(isoreader)
library(QCIRMS)
dxf_file <- "170506_NaHCO3 L + NaCl U_.dxf"
dataPath <- system.file("extdata/dxf_files/abiotic", dxf_file, package = "QCIRMS")
file.summ<-read_summary(dataPath)
```
### Abstract

#### Contact
[lily-clough@utulsa.edu](lily-clough@utulsa.edu)
[brett-mckinney@utulsa.edu](brett-mckinney@utulsa.edu)
