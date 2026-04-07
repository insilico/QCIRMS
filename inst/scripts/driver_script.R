#pak::pak("insilico/QCIRMS")
library(QCIRMS)

setwd("~/Documents/Papers/qcirms/test")

out <- render_qcirms_report(
  # make sure dirs don't have spaces
  DATA_PATH = "220318_1pct_CO2_NaCl_NaSO4_brines",
  REF_TIME_PATH = "220318_1pct_CO2_NaCl_NaSO4_brines/qc",
  check_internal_standards = T, 
  PLOTS_PATH = "test_plots",
  RESULT_PATH = "test_results",
  output_format = "html", # pdf
  report_rmd = "test_report.Rmd",
  report_title = "Test Report"
)


