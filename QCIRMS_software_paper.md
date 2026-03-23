# QCIRMS: an R workflow for quality control of continuous-flow IRMS DXF data

**Brett McKinney** and **Lily Clough**  
The University of Tulsa

## Abstract

Continuous-flow isotope ratio mass spectrometry (IRMS) experiments often generate large
numbers of vendor DXF files whose downstream quality control is partly manual, partly
ad hoc, and difficult to reproduce across experiments. QCIRMS is an R package that
provides a consistent workflow for reading DXF files, extracting vendor tables and raw
chromatograms, plotting full experiment batches, and applying directory-level quality
control filters to reference peaks, sample peaks, and internal standards. The package
builds on `isoreader` for file parsing and adds IRMS-specific utilities for expected
reference-peak timing checks, isotopic precision checks, relative intensity checks, and
export of curated outputs for downstream analysis. The result is a transparent and
scriptable QA/QC pipeline that reduces manual review time and improves reproducibility
for stable-isotope workflows.

## Introduction

Stable-isotope workflows are highly sensitive to experimental drift, peak misassignment,
retention-time variation, and inconsistent handling of standards. In many laboratories,
quality control still relies on a mixture of instrument software, spreadsheet review,
and investigator-specific conventions. This slows analysis and makes it harder to apply
the same QC rules across studies or operators. An R-based workflow is attractive because
it can unify file import, visualization, filtering, and result export in a single,
version-controlled environment.

QCIRMS was developed to address this need for continuous-flow IRMS data stored in Thermo
DXF format. Rather than replacing the vendor export, QCIRMS wraps it in a reproducible
analysis workflow. The package supports inspection of individual files, batch processing
of experiment directories, and an end-to-end wrapper, `QAQC_IRMS()`, for quality control
of complete IRMS runs.

## Software overview

QCIRMS provides two complementary levels of functionality. The first level is exploratory.
Functions such as `read_summary()`, `vendor_info()`, `raw_data()`, `reference_values_ratio()`,
and `reference_values_no_ratio()` allow users to inspect a single DXF file and recover
the key vendor tables needed for method development, troubleshooting, and teaching. Plotting
utilities such as `generic_raw_plot()` and `generic_plot_all_raw()` enable rapid visual
review of chromatograms for one file or for an entire directory.

The second level is analytical and workflow-oriented. QCIRMS includes utilities to list
DXF files, read directories into batch objects, compare observed reference-peak times to
expected values, and screen analyses based on isotopic variability and peak-intensity
consistency. The main wrapper, `QAQC_IRMS()`, applies these steps to a directory of files
and returns curated reference and sample peak tables, with optional internal-standard
calibration output. This structure makes the package useful both for exploratory laboratory
review and for reproducible post-acquisition processing.

## Example workflow

A typical workflow begins by reading a representative DXF file to inspect the vendor
table and raw chromatogram. Users can then point QCIRMS to a directory of DXF files,
generate a batch PDF of raw chromatograms, and verify the expected times of reference
peaks from a reference table. Once these inputs are available, `QAQC_IRMS()` performs
directory-level QA/QC using user-defined thresholds for the standard deviations of
reference and sample isotope ratios, expected non-sample peak counts, and relative
differences in peak intensity. The package writes tabular outputs for downstream analysis
while retaining the intermediate information needed to audit the filtering decisions.

This design is especially useful in studies with many replicate runs or multiple brine
chemistries, where manual review is time-consuming and difficult to standardize. Because
the workflow is scripted, it can be re-run when thresholds change, new data are added,
or package code is updated.

## Availability and impact

QCIRMS is distributed as an R package through GitHub and currently depends on the GitHub
package `isoreader` for DXF parsing. The software is intended for researchers who need a
lightweight, transparent workflow for IRMS QA/QC in R, especially when they want direct
access to raw chromatograms and vendor-derived peak tables rather than only final isotope
summaries. By combining file import, visualization, rule-based QC, and export in one
environment, QCIRMS lowers the friction of reproducible IRMS data processing and provides
a foundation for future extensions such as vignettes, automated reporting, and improved
support for internal-standard calibration.

## Software availability

- Name: QCIRMS
- Language: R
- Source code: GitHub repository for QCIRMS
- Dependency for DXF parsing: `isoreader`
- Suggested tutorial format: R Markdown / HTML vignette
