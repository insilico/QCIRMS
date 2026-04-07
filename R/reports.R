#' Write and render a QCIRMS report
#'
#' Convenience wrapper around [write_qcirms_report_rmd()] and
#' [rmarkdown::render()]. This function writes a parameterized QCIRMS report
#' template and immediately renders it to HTML or PDF.
#'
#' @param DATA_PATH Character string giving the directory containing input
#'   `.dxf` files.
#' @param REF_TIME_PATH Character string giving the directory containing the
#'   expected reference-times file.
#' @param PLOTS_PATH Character string giving the directory where plot files
#'   created by the report should be written.
#' @param RESULT_PATH Character string giving the directory where result files,
#'   the generated `.Rmd`, and the rendered report should be written.
#' @param output_format Character string giving the render format. Supported
#'   values are `"html_document"` and `"pdf_document"`.
#' @param output_file Optional character string giving the rendered report
#'   filename. If `NULL`, a default is chosen based on `output_format`.
#' @param report_rmd Optional character string giving the `.Rmd` filename to
#'   write. Defaults to `file.path(RESULT_PATH, "QCIRMS_report.Rmd")`.
#' @param report_title Character string giving the report title.
#' @param ref_file Character string giving the expected reference-times filename.
#' @param check_internal_standards Logical; whether the report should run
#'   `QAQC_IRMS()` with internal-standard checking enabled.
#' @param overwrite Logical; whether to overwrite an existing `.Rmd` file.
#' @param clean Logical; passed to [rmarkdown::render()]. If `TRUE`, intermediate
#'   files produced during rendering are cleaned up.
#' @param envir Environment in which to render the report. Defaults to
#'   `new.env(parent = globalenv())`.
#'
#' @return Invisibly returns a list with elements:
#' \itemize{
#'   \item `rmd_file`: path to the written `.Rmd`
#'   \item `rendered_file`: path to the rendered output
#' }
#'
#' @examples
#' \dontrun{
#' out <- render_qcirms_report(
#'   DATA_PATH = "inst/extdata/dxf_files/abiotic",
#'   REF_TIME_PATH = "inst/extdata/qc",
#'   PLOTS_PATH = "plots",
#'   RESULT_PATH = "results",
#'   output_format = "html_document"
#' )
#'
#' out$rendered_file
#' }
#'
#' @export
render_qcirms_report <- function(
    DATA_PATH,
    REF_TIME_PATH,
    PLOTS_PATH,
    RESULT_PATH,
    output_format = c("html_document", "pdf_document"),
    output_file = NULL,
    report_rmd = file.path(RESULT_PATH, "QCIRMS_report.Rmd"),
    report_title = "QCIRMS Processing Report",
    ref_file = "referencePeaks_expectedTimes.txt",
    check_internal_standards = TRUE,
    overwrite = TRUE,
    clean = TRUE,
    envir = new.env(parent = globalenv())
) {
  output_format <- match.arg(output_format)
  
  dir.create(PLOTS_PATH, recursive = TRUE, showWarnings = FALSE)
  dir.create(RESULT_PATH, recursive = TRUE, showWarnings = FALSE)
  
  yaml_formats <- switch(
    output_format,
    html_document = "html",
    pdf_document = "pdf"
  )
  
  rmd_file <- write_qcirms_report_rmd(
    DATA_PATH = DATA_PATH,
    REF_TIME_PATH = REF_TIME_PATH,
    PLOTS_PATH = PLOTS_PATH,
    RESULT_PATH = RESULT_PATH,
    output_rmd = report_rmd,
    report_title = report_title,
    ref_file = ref_file,
    check_internal_standards = check_internal_standards,
    output_formats = yaml_formats,
    overwrite = overwrite
  )
  
  if (is.null(output_file)) {
    output_file <- switch(
      output_format,
      html_document = "QCIRMS_report.html",
      pdf_document = "QCIRMS_report.pdf"
    )
  }
  
  rendered_file <- rmarkdown::render(
    input = rmd_file,
    output_format = output_format,
    output_file = output_file,
    output_dir = RESULT_PATH,
    envir = envir,
    clean = clean
  )
  
  invisible(list(
    rmd_file = rmd_file,
    rendered_file = rendered_file
  ))
}

#' Write a QCIRMS report template as an R Markdown file
#'
#' Creates an `.Rmd` report template based on the QCIRMS batch-processing
#' workflow. The generated report can be knitted to HTML or PDF and is designed
#' to both write output files to disk and embed key tables, summaries, and plots
#' directly in the report.
#'
#' The report includes:
#' \itemize{
#'   \item path summary
#'   \item DXF file discovery
#'   \item single-file preview
#'   \item reference-time table
#'   \item batch raw-data summary
#'   \item batch plotting
#'   \item full `QAQC_IRMS()` execution
#'   \item QA/QC output summaries
#'   \item optional internal-standard summaries
#'   \item session information
#' }
#'
#' @param DATA_PATH Character string giving the directory containing input
#'   `.dxf` files.
#' @param REF_TIME_PATH Character string giving the directory containing the
#'   expected reference-times file.
#' @param PLOTS_PATH Character string giving the directory where plot files
#'   created by the report should be written.
#' @param RESULT_PATH Character string giving the directory where result files
#'   and the generated `.Rmd` should be written.
#' @param output_rmd Character string giving the path of the output `.Rmd` file.
#'   Defaults to `file.path(RESULT_PATH, "QCIRMS_report.Rmd")`.
#' @param report_title Character string giving the report title.
#' @param ref_file Character string giving the expected reference-times filename.
#'   Defaults to `"referencePeaks_expectedTimes.txt"`.
#' @param check_internal_standards Logical; whether the generated report should
#'   run `QAQC_IRMS()` with internal-standard checking enabled.
#' @param output_formats Character vector of output formats to include in the
#'   YAML header. Supported values are `"html"` and `"pdf"`.
#' @param overwrite Logical; whether to overwrite an existing `.Rmd` file.
#'
#' @return Invisibly returns the path to the generated `.Rmd` file.
#'
#' @details
#' The generated report is parameterized through YAML `params`, but the file is
#' written with the provided path values already filled in. Users can knit it
#' directly in RStudio or render it programmatically with `rmarkdown::render()`.
#'
#' The report writes results to `PLOTS_PATH` and `RESULT_PATH`, while also
#' embedding selected previews and summary tables in the knitted document.
#'
#' @examples
#' \dontrun{
#' report_file <- write_qcirms_report_rmd(
#'   DATA_PATH = "inst/extdata/dxf_files/abiotic",
#'   REF_TIME_PATH = "inst/extdata/qc",
#'   PLOTS_PATH = "plots",
#'   RESULT_PATH = "results"
#' )
#'
#' rmarkdown::render(report_file, output_format = "html_document")
#' }
#'
#' @export
write_qcirms_report_rmd <- function(
    DATA_PATH,
    REF_TIME_PATH,
    PLOTS_PATH,
    RESULT_PATH,
    output_rmd = file.path(RESULT_PATH, "QCIRMS_report.Rmd"),
    report_title = "QCIRMS Processing Report",
    ref_file = "referencePeaks_expectedTimes.txt",
    check_internal_standards = TRUE,
    output_formats = c("html", "pdf"),
    overwrite = TRUE
) {
  stopifnot(is.character(DATA_PATH), length(DATA_PATH) == 1L)
  stopifnot(is.character(REF_TIME_PATH), length(REF_TIME_PATH) == 1L)
  stopifnot(is.character(PLOTS_PATH), length(PLOTS_PATH) == 1L)
  stopifnot(is.character(RESULT_PATH), length(RESULT_PATH) == 1L)
  stopifnot(is.character(output_rmd), length(output_rmd) == 1L)
  stopifnot(is.character(report_title), length(report_title) == 1L)
  stopifnot(is.character(ref_file), length(ref_file) == 1L)
  stopifnot(is.logical(check_internal_standards), length(check_internal_standards) == 1L)
  stopifnot(is.logical(overwrite), length(overwrite) == 1L)
  
  output_formats <- unique(tolower(output_formats))
  allowed_formats <- c("html", "pdf")
  bad_formats <- setdiff(output_formats, allowed_formats)
  if (length(bad_formats) > 0L) {
    stop("Unsupported output_formats: ", paste(bad_formats, collapse = ", "))
  }
  if (length(output_formats) == 0L) {
    stop("At least one output format must be requested.")
  }
  
  if (!dir.exists(DATA_PATH)) {
    stop("DATA_PATH does not exist: ", DATA_PATH)
  }
  if (!dir.exists(REF_TIME_PATH)) {
    stop("REF_TIME_PATH does not exist: ", REF_TIME_PATH)
  }
  
  dir.create(PLOTS_PATH, recursive = TRUE, showWarnings = FALSE)
  dir.create(RESULT_PATH, recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(output_rmd) && !overwrite) {
    stop("output_rmd already exists and overwrite = FALSE: ", output_rmd)
  }
  
  ref_full_path <- file.path(REF_TIME_PATH, ref_file)
  if (!file.exists(ref_full_path)) {
    warning("Reference-times file does not currently exist: ", ref_full_path)
  }
  
  normalize_slash <- function(x) {
    gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
  }
  
  yaml_lines <- c(
    "---",
    paste0('title: "', report_title, '"'),
    'author: "QCIRMS"',
    'date: "`r Sys.Date()`"',
    "output:"
  )
  
  if ("html" %in% output_formats) {
    yaml_lines <- c(
      yaml_lines,
      "  html_document:",
      "    toc: true",
      "    toc_depth: 3",
      "    df_print: paged"
    )
  }
  
  if ("pdf" %in% output_formats) {
    yaml_lines <- c(
      yaml_lines,
      "  pdf_document:",
      "    toc: true",
      "    number_sections: true"
    )
  }
  
  yaml_lines <- c(
    yaml_lines,
    "params:",
    paste0('  DATA_PATH: "', normalize_slash(DATA_PATH), '"'),
    paste0('  REF_TIME_PATH: "', normalize_slash(REF_TIME_PATH), '"'),
    paste0('  PLOTS_PATH: "', normalize_slash(PLOTS_PATH), '"'),
    paste0('  RESULT_PATH: "', normalize_slash(RESULT_PATH), '"'),
    paste0('  REF_FILE: "', ref_file, '"'),
    paste0("  CHECK_INTERNAL_STANDARDS: ", if (check_internal_standards) "true" else "false"),
    "always_allow_html: true",
    "---",
    ""
  )
  
  results_prefix <- basename(normalize_slash(DATA_PATH))
  
  body_lines <- c(
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(",
    "  echo = TRUE,",
    "  warning = FALSE,",
    "  message = FALSE,",
    "  fig.align = 'center'",
    ")",
    "library(QCIRMS)",
    "library(isoreader)",
    "library(dplyr)",
    "library(knitr)",
    "```",
    "",
    "# Overview",
    "",
    "This report was generated automatically from a QCIRMS processing template.",
    "It follows the batch-processing structure of the QCIRMS workflow and includes both written output files and embedded summaries.",
    "",
    "# Paths",
    "",
    "```{r paths}",
    "DATA_PATH <- params$DATA_PATH",
    "REF_TIME_PATH <- params$REF_TIME_PATH",
    "PLOTS_PATH <- params$PLOTS_PATH",
    "RESULT_PATH <- params$RESULT_PATH",
    "REF_FILE <- params$REF_FILE",
    "CHECK_INTERNAL_STANDARDS <- params$CHECK_INTERNAL_STANDARDS",
    "",
    "dir.create(PLOTS_PATH, recursive = TRUE, showWarnings = FALSE)",
    "dir.create(RESULT_PATH, recursive = TRUE, showWarnings = FALSE)",
    "",
    "path_tbl <- data.frame(",
    "  Setting = c('DATA_PATH', 'REF_TIME_PATH', 'PLOTS_PATH', 'RESULT_PATH', 'REF_FILE', 'CHECK_INTERNAL_STANDARDS'),",
    "  Value = c(DATA_PATH, REF_TIME_PATH, PLOTS_PATH, RESULT_PATH, REF_FILE, CHECK_INTERNAL_STANDARDS),",
    "  stringsAsFactors = FALSE",
    ")",
    "knitr::kable(path_tbl, caption = 'Report parameters')",
    "```",
    "",
    "# Discover DXF files",
    "",
    "```{r discover-files}",
    "fileNames <- all_dxf_files(path = DATA_PATH)",
    "cat('Number of DXF files found:', length(fileNames), '\\n')",
    "if (length(fileNames) > 0) {",
    "  knitr::kable(head(data.frame(file = fileNames, stringsAsFactors = FALSE), 10),",
    "               caption = 'First DXF files discovered')",
    "} else {",
    "  stop('No .dxf files found in DATA_PATH.')",
    "}",
    "```",
    "",
    "# Single-file preview",
    "",
    "```{r single-file-preview}",
    "first_file <- fileNames[1]",
    "cat('Preview file:', first_file, '\\n')",
    "",
    "vend.df <- vendor_info(file.path(DATA_PATH, first_file))",
    "raw.df  <- raw_data(file.path(DATA_PATH, first_file))",
    "",
    "knitr::kable(head(vend.df[, seq_len(min(10, ncol(vend.df))), drop = FALSE]),",
    "             caption = 'Vendor information preview')",
    "knitr::kable(head(raw.df), caption = 'Raw chromatogram preview')",
    "",
    "write.csv(vend.df, file.path(RESULT_PATH, 'single_file_vendor_info.csv'), row.names = FALSE)",
    "write.csv(raw.df,  file.path(RESULT_PATH, 'single_file_raw_data.csv'), row.names = FALSE)",
    "```",
    "",
    "```{r single-file-plot, fig.width=8, fig.height=5}",
    "generic_raw_plot(raw.df = raw.df, title = first_file)",
    "```",
    "",
    "# Reference information",
    "",
    "```{r reference-info}",
    "expRef.df <- read.table(file.path(REF_TIME_PATH, REF_FILE), header = TRUE)",
    "knitr::kable(expRef.df, caption = 'Expected reference peak times')",
    "write.csv(expRef.df, file.path(RESULT_PATH, 'expected_reference_times.csv'), row.names = FALSE)",
    "```",
    "",
    "# Batch raw-data processing",
    "",
    "```{r batch-raw}",
    "rawList <- raw_data_all(fileNames, path = DATA_PATH)",
    "cat('Number of raw chromatogram objects:', length(rawList), '\\n')",
    "raw_lengths <- vapply(rawList, function(x) if (is.null(x)) 0L else nrow(x), integer(1))",
    "raw_tbl <- data.frame(file = fileNames, n_rows = raw_lengths, stringsAsFactors = FALSE)",
    "knitr::kable(head(raw_tbl, 10), caption = 'Raw data row counts')",
    "write.csv(raw_tbl, file.path(RESULT_PATH, 'raw_data_row_counts.csv'), row.names = FALSE)",
    "```",
    "",
    "# Batch plots",
    "",
    "```{r batch-plots}",
    "generic_plot_all_raw(",
    "  raw.list = rawList,",
    "  path = PLOTS_PATH,",
    "  write_pdf = TRUE,",
    paste0("  pdf_name = '", results_prefix, "_all_dxf.pdf'"),
    ")",
    paste0("cat('Batch plot PDF written to:', file.path(PLOTS_PATH, '", results_prefix, "_all_dxf.pdf'), '\\n')"),
    "```",
    "",
    "# Run QCIRMS QA/QC",
    "",
    "```{r qaqc-run}",
    "start_time <- Sys.time()",
    "currFiltered <- QAQC_IRMS(",
    "  unfilteredPath = DATA_PATH,",
    "  expRef.df = expRef.df,",
    "  checkIntStand = CHECK_INTERNAL_STANDARDS,",
    "  internalStandID = c('L1', 'H1', 'LW'),",
    paste0("  dataName = '", results_prefix, "',"),
    "  maxPkNum = 18,",
    "  expectedNonSampPks = 7,",
    "  sdCrefIso.thresh = 0.1,",
    "  sdOrefIso.thresh = 0.1,",
    "  checkRelDiffIntensity = TRUE,",
    "  refSamp_relDiffInt.thresh = 0.75,",
    "  amplName = 'Ampl44',",
    "  ref_relDiffInt.thresh = 0.1,",
    "  sdCsampIso.thresh = 0.3,",
    "  sdOsampIso.thresh = 0.2,",
    "  outPath = RESULT_PATH,",
    "  verbose = TRUE",
    ")",
    "elapsed <- Sys.time() - start_time",
    "cat('QA/QC runtime:', elapsed, '\\n')",
    "```",
    "",
    "# QA/QC outputs",
    "",
    "```{r qaqc-outputs}",
    "refs.dat <- currFiltered[[1]]",
    "samps.dat <- currFiltered[[2]]",
    "",
    "cat('Reference rows after QA/QC:', nrow(refs.dat), '\\n')",
    "cat('Sample rows after QA/QC:', nrow(samps.dat), '\\n')",
    "",
    "knitr::kable(head(refs.dat), caption = 'Reference peaks passing QA/QC')",
    "knitr::kable(head(samps.dat), caption = 'Sample peaks passing QA/QC')",
    "",
    "write.csv(refs.dat, file.path(RESULT_PATH, 'report_refs_passed.csv'), row.names = FALSE)",
    "write.csv(samps.dat, file.path(RESULT_PATH, 'report_samples_passed.csv'), row.names = FALSE)",
    "",
    "summary_tbl <- data.frame(",
    "  object = c('reference_peaks', 'sample_peaks'),",
    "  n_rows = c(nrow(refs.dat), nrow(samps.dat)),",
    "  n_cols = c(ncol(refs.dat), ncol(samps.dat)),",
    "  stringsAsFactors = FALSE",
    ")",
    "knitr::kable(summary_tbl, caption = 'Summary of QA/QC outputs')",
    "write.csv(summary_tbl, file.path(RESULT_PATH, 'report_summary.csv'), row.names = FALSE)",
    "```",
    "",
    "# Internal standards",
    "",
    "```{r internal-standards}",
    "if (CHECK_INTERNAL_STANDARDS && length(currFiltered) >= 3) {",
    "  int_stand.list <- currFiltered[[3]]",
    "",
    "  if (length(int_stand.list) >= 2 && !is.null(int_stand.list[[2]])) {",
    "    knitr::kable(int_stand.list[[2]], caption = 'Average internal-standard deltas')",
    "    write.csv(int_stand.list[[2]], file.path(RESULT_PATH, 'internal_standard_means.csv'), row.names = FALSE)",
    "  }",
    "",
    "  if (length(int_stand.list) >= 3 && !is.null(int_stand.list[[3]])) {",
    "    knitr::kable(int_stand.list[[3]], caption = 'Internal-standard standard deviations')",
    "    write.csv(int_stand.list[[3]], file.path(RESULT_PATH, 'internal_standard_sds.csv'), row.names = FALSE)",
    "  }",
    "} else {",
    "  cat('Internal standards were not requested or were not returned.\\n')",
    "}",
    "```",
    "",
    "# Session information",
    "",
    "```{r session-info}",
    "sessionInfo()",
    "```",
    ""
  )
  
  writeLines(c(yaml_lines, body_lines), con = output_rmd)
  invisible(output_rmd)
}

