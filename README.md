# QCIRMS
Quality Control Analysis for Isotope Ratio Mass Spectrometry Data

### To install:

```
#install.packages("pak")
library(pak)
pak::pak("insilico/QCIRMS")
```
If the above doesn't work, try this.
```
devtools::install_github("isoverse/isoreader")
devtools::install_github("insilico/QCIRMS")
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

This is the QCIRMS package, which is designed to perform quality control 
analysis on isotope ratio mass spectrometry (IRMS) data. 
The package provides tools for reading, processing, and analyzing IRMS 
data files, allowing researchers to assess the quality of their measurements 
and identify potential issues.

#### Contact
[lily-clough@utulsa.edu](lily-clough@utulsa.edu)
[brett-mckinney@utulsa.edu](brett-mckinney@utulsa.edu)
