

temp_data_frame <- data.frame(
    fileID = "test_dataframe_file.dxf",
    Identifier1 = "sample",
    Analysis = 1,
    Preparation = "prep",
    DataTime = "date",
    PeakNr = 1:8,
    Nr. = 1:8,
    Start = c(10, 20, 30, 40, 50, 60, 70, 80),
    Rt = c(11, 21, 31, 41, 51, 61, 71, 81),
    End = c(12, 22, 32, 42, 52, 62, 72, 82),
    Ampl44 = c(100, 101, 100, 101, 100, 100, 101, 100),
    Ampl45 = 1,
    Ampl46 = 1,
    d13C12C = c(rep(1, 5), rep(2, 3)),
    d18O16O = c(rep(1, 5), rep(2, 3)),
    IsRef = c(rep(TRUE, 5), rep(FALSE, 3)),
    stringsAsFactors = FALSE
  )

# 
# make_exp_ref_df <- function(vend.df = make_vend_df()) {
#   data.frame(
#     Ref_Peak_Nr = 1:8,
#     Expected_Start = c(vend.df$Start[1:5], 999, 999, 999),
#     Expected_Rt = c(vend.df$Rt[1:5], 999, 999, 999),
#     Expected_End = c(vend.df$End[1:5], 999, 999, 999)
#   )
# }
# 
# make_ref_times_output <- function(vend.df = make_vend_df(), expRef.df = make_exp_ref_df(vend.df)) {
#   QCIRMS:::reference_times_checkDXF(
#     vend.df = vend.df,
#     expectedPeak.num = 5,
#     diff.t = 1,
#     expectedStart = expRef.df$Expected_Start,
#     expectedRt = expRef.df$Expected_Rt,
#     expectedEnd = expRef.df$Expected_End,
#     verbose = FALSE
#   )
# }
# 
# make_processed_list <- function() {
#   vend.df <- make_vend_df()
#   expRef.df <- make_exp_ref_df(vend.df)
#   QCIRMS:::removeFailedAnalysesDXF(
#     sepList = list(vend.df),
#     expRef.df = expRef.df,
#     maxPkNum = 10,
#     expRefPkNr = 5,
#     expectedNonSampPks = 5,
#     diff.t = 1,
#     sdCrefIso.thresh = 0.1,
#     sdOrefIso.thresh = 0.1,
#     relDiffInt.thresh = 0.1,
#     amplName = "Ampl44",
#     sdCsampIso.thresh = 0.1,
#     sdOsampIso.thresh = 0.1,
#     filter_flushPk = FALSE,
#     filter_first_sampPk = FALSE,
#     flushExpT = 60,
#     flushTint = 5,
#     firstSampExpT = 70,
#     firstSampTint = 5,
#     firstSampExpPkNr = 6,
#     verbose = FALSE
#   )
# }

test_that("function returns expected data type", {
  result <- separate_by_analysis_numDXF(temp_data_frame)

  expect_type(result, "list")
})
