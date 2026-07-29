suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(isoreader))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pracma))


#file settings

#removed system.file, since it was directing to the files in the QCIRMS library instead of my personal copy

#testfiles_dir: setwd("/Users/abigailblakey/Documents/git/QCIRMS_my_copy_july")

#input for single builtin .dxf file, test_path allows tests to access file path

builtin_file <- normalizePath(testthat::test_path(
  "fixtures",
  "test_170506_NaHCO3 L + NaCl U.dxf"
), mustWork = TRUE)

#input for dxf dataframe (vend.df and raw_data.df)

suppressMessages(invisible(capture.output(vend_builtin_file <- vendor_info(builtin_file))))
suppressMessages(invisible(capture.output(raw_builtin_file <- raw_data(builtin_file))))

#input for exp.df

ref_time_file <- normalizePath(testthat::test_path(
  "fixtures",
  "referencePeaks_expectedTimes copy.txt"
), mustWork = TRUE)

test_exp_df <- read.table(ref_time_file)

#input for folder of .dxf files

builtin_dxf_folder <- normalizePath(testthat::test_path(
  "fixtures",
  "test_dxf_folder"
), mustWork = TRUE)

#vector list of multiple .dxf files

builtin_dxf_folder_vector_list <- list.files(
  path = builtin_dxf_folder,
  pattern = "\\.dxf$",
  full.names = TRUE,
  ignore.case = TRUE
)

#input for .dxf folder vendor/raw data

suppressMessages(invisible(capture.output(vend_builtin_dxf_folder <- vendor_info_all(builtin_dxf_folder))))
suppressMessages(invisible(capture.output(raw_builtin_dxf_folder <- raw_data_all(builtin_dxf_folder_vector_list))))

#analysis numbers corresponding to the 4 .dxf files in builtin_dxf_folder
test_analysisNum.vec <- c(11890,11869,11891,11911)

#full rows from checkSamps for the 4 .dxf files in builtin_dxf_folder

check_file <- normalizePath(testthat::test_path(
  "fixtures",
  "220318_1pct_CO2_NaCl_NaSO4_brines_samps_allChecks copy.csv"
), mustWork = TRUE)

check_df <- read.csv(check_file)

filtered_check_df <- check_df %>%
  dplyr::filter(fileID %in% c(
    "0_1 m Na2SO4_(1).dxf", #fails
    "6_0 m NaCl_(2).dxf", #passes
    "5_0 m NaCl_.dxf", #passes
    "2_0 m Na2SO4_.dxf" #passes
  ))
  

#input for check file (plot_pass_fail_spectra())





# vector of analysis numbers

#processed list

#processed.list <- removeFailedAnalysesDXF(sepList=combList,
#'                                           expRef.df=expRef.df,
#'                                           maxPkNum=18,
#'                                           expRefPkNr=dim(expRef.df)[1],
#'                                           diff.t=10,
#'                                           sdCrefIso.thresh=0.1, 
#'                                           expectedNonSampPks=7,
#'                                           relDiffInt.thresh=0.1, 
#'                                           amplName="Ampl44",
#'                                           sdOrefIso.thresh=0.1,
#'                                           sdCsampIso.thresh=0.3,
#'                                           sdOsampIso.thresh=0.2,
#'                                           flushExpT=135,
#'                                           flushTint=15,
#'                                           firstSampExpT=275,
#'                                           firstSampTint=15,
#'                                           verbose=T)# temp_brine_test_folder_path <- testthat::test_path(
#   "..",
#   "..",
#   "220318_1pct_CO2_NaCl_NaSO4_brines"
# )

#brine_test_folder <- "220318_1pct_CO2_NaCl_NaSO4_brines"

#check if all files exist

# result <- expect_true(file.exists(builtin_file))
# 
# if (result == TRUE) {
#   print("file exists!")
# }
