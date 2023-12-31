% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QAQC_IRMS_functions.R
\name{find_files_by_analysis_num}
\alias{find_files_by_analysis_num}
\title{find_files_by_analysis_num}
\usage{
find_files_by_analysis_num(path, analysisNum.vec)
}
\arguments{
\item{path}{to dxf files to search}

\item{analysisNum.vec}{vector of analysis numbers}
}
\value{
vector containing file names of dxf files for the specified analysis numbers
}
\description{
find_files_by_analysis_num
}
\examples{
Usage example
outlier_analysis.vec <- outlier_test_meta.dat$Analysis 
abiotic_path <- "./data/dxf_files/abiotic/"
abio_outlier_dxf_files <- find_files_by_analysis_num(path = abiotic_path, analysisNum.vec = outlier_analysis.vec)
}
