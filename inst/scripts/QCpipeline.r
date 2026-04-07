#devtools::install_github("isoverse/isoreader")
library(isoreader)
library(dplyr)
library(pracma)

# set working directory - sets wd of this script wherever you saved it 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load functions
#source("QAQC_IRMS_functions_v2.R")
library(QCIRMS)

##################

### basics - working with one dxf file 

# Can print summaries of any of the data files - info on contents

# example file in directory for library from Github (file_path)
dxf_file <- "170506_NaHCO3 L + NaCl U_.dxf"
file_path <- system.file("extdata/dxf_files/abiotic", dxf_file, package = "QCIRMS")
file.summ<-read_summary(file_path)
#                   Length Class           Mode
# version            1     package_version list
# read_options       4     -none-          list
# file_info         16     tbl_df          list
# method_info        3     -none-          list
# raw_data           5     tbl_df          list
# vendor_data_table 39     tbl_df          list


# Can get isoreader file information - labels input at the time of the experiment
(file_info.df<-file_info(files=file_path))
#                        file_id       Identifier1 Analysis             Preparation       Date_and_Time
#1 170506_NaHCO3 L + NaCl U_.dxf NaHCO3 L + NaCl U     3355 0.115ml CO2 in He 24hrs 2017-05-06 15:40:50


# Can get isoreader vendor info - "mass spectrometry" data - IRMS measurements and derived values
vend.df<-vendor_info(file_path)
head(vend.df)[,1:10]
#                        file_id Nr.   Start      Rt     End  Ampl 44  Ampl 45  Ampl 46   BGD 44    BGD 45
#1 170506_NaHCO3 L + NaCl U_.dxf   1  27.170  47.443  50.787 2425.983 2788.829 3307.896 1.207758 0.5735142
#2 170506_NaHCO3 L + NaCl U_.dxf   2  66.880  87.153  90.497 2425.019 2787.828 3306.349 1.599683 0.9146100
#3 170506_NaHCO3 L + NaCl U_.dxf   3 131.043 133.551 140.030 3592.126 4210.890 4937.000 1.429455 0.7790426
dim(vend.df)
# [1] 17 40 # not all have data 


# Can get the raw detection data (intensity vs time chromatogram)
raw.df<-raw_data(file_path)
head(raw.df)
#                         file_id tp time.s   v44.mV    v45.mV   v46.mV
# 1 170506_NaHCO3 L + NaCl U_.dxf  1  0.209 1.128612 0.4017511 2.038189
# 2 170506_NaHCO3 L + NaCl U_.dxf  2  0.418 1.130524 0.3998433 1.996066
# 3 170506_NaHCO3 L + NaCl U_.dxf  3  0.627 1.117143 0.4208301 1.942465
# 4 170506_NaHCO3 L + NaCl U_.dxf  4  0.836 1.117143 0.3941199 1.942465
# 5 170506_NaHCO3 L + NaCl U_.dxf  5  1.045 1.124789 0.4151062 1.927152
# 6 170506_NaHCO3 L + NaCl U_.dxf  6  1.254 1.113320 0.4131983 1.967350
dim(raw.df)
# [1] 4298    6


# Can get the resistor information - sometimes useful according to Bethany
(resist.df<-resistor_data(file_path))
#                         file_id cup R.Ohm mass
# 1 170506_NaHCO3 L + NaCl U_.dxf   1 3e+08   44
# 2 170506_NaHCO3 L + NaCl U_.dxf   2 3e+10   45
# 3 170506_NaHCO3 L + NaCl U_.dxf   3 1e+11   46


# Can get isotopic reference values with ratios ...
(stand_ratio.df<-reference_values_ratio(file_path))
#                         file_id standard gas delta_name delta_value reference element ratio_name ratio_value
# 1 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 13C/12C       -36.9      VPDB       C  R 13C/12C  0.01118020
# 2 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 13C/12C       -36.9      VPDB       O  R 18O/16O  0.00206720
# 3 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 13C/12C       -36.9      VPDB       O  R 17O/16O  0.00038600
# 4 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 18O/16O       -40.0     VSMOW       H    R 2H/1H  0.00015575
# 5 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 18O/16O       -40.0     VSMOW       O  R 17O/16O  0.00037990
# 6 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 18O/16O       -40.0     VSMOW       O  R 18O/16O  0.00200520

# ... and without ratios
(stand_no_ratio.df<-reference_values_no_ratio(file_path))
#                         file_id standard gas delta_name delta_value reference
# 1 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 13C/12C       -36.9      VPDB
# 2 170506_NaHCO3 L + NaCl U_.dxf CO2_zero CO2  d 18O/16O       -40.0     VSMOW


# Can also get peak areas
(trap_area_allPks(raw.df,vend.df,mV.rawName="v44.mV"))
#    Pk_Nr  trap_area
# 1      1 47.6436699
# 2      2 47.6571723
# 3      3  5.6953350
# 4      4 47.6760681
# 5      5 47.4243105
# 6      6 23.7051077
# 7      7  0.2369494
# 8      8 22.5549720
# 9      9 21.4485258
# 10    10 20.3915933
# 11    11 19.3905900
# 12    12 18.4485201
# 13    13 17.5371007
# 14    14 16.6757975
# 15    15 15.7937530
# 16    16 15.0370583
# 17    17 47.4650624


# Can plot raw data (Intensity (mV) vs Rt)
generic_raw_plot(raw.df=raw.df,title=dxf_file)

## save plot to pdf in a results directory with the data label
plotsPath<-"./QCresults/abiotic/" # could make a directory for results and store plots there
generic_raw_plot(raw.df=raw.df, # intensity vs time dataframe
                 title=dxf_file, # title on plot
                 path=plotsPath, # path to write file
                 write_pdf=T, # write pdf
                 pdf_name="one_dxf.pdf") 



# Can plot raw data of all files in a vector or directory
# first get all filenames of .dxf files in the data directory

# path to dxf files (raw IRMS experimental data)
dataPath<-"./data/dxf_files/abiotic/" # could save your data in a structure like this
# get dxf file names
fileNames<-all_dxf_files(path=dataPath) # looks in current wd if path not specified 
# make list of raw ampl vs time dataframes (time-series)
rawList<-raw_data_all(fileNames,path=dataPath) # looks in curr wd if path to data not specified 
length(rawList)
# [1] 789

# plot all raw data using code that writes pdfs in the specified working directory
generic_plot_all_raw(rawList, # list of dataframes
                     path=plotsPath, # use path for plots defined before
                     write_pdf=T,
                     pdf_name="all_dxf.pdf") 




###################
### QA/QC for a directory of volatile CO2 IRMS experiments

# get expected times for reference peaks
# path to ref peaks data in library, or set to where you stored this file
refFile<-"referencePeaks_expectedTimes.txt"
refTime_path <- system.file("extdata/qc/", refFile, package = "QCIRMS")

expRef.df<-read.table(refTime_path)
expRef.df
#   Ref_Peak_Nr Expected_Start Expected_Rt Expected_End
# 1           1        27.1700     47.3888      50.4929
# 2           2        67.0193     87.1840      90.3422
# 3           4       166.5653    186.7608     189.8572
# 4           5       206.4920    226.5096     229.6987
# 5          16       843.3227    863.5029     866.6069


# provide dataset name
dataName<-"abiotic"

# paths to data and results... make sure slash is part of results path to write files
# same directory plots are in, now store QC results
resultsPath<-"./QCresults/abiotic/"

# run QA/QC on directory of dxf files 
# execute start <- ... to QAtime together (select then run) to print runtime
start<-Sys.time()
currFiltered<-QAQC_IRMS(unfilteredPath=dataPath, # specified previously
                        expRef.df=expRef.df, 
                        checkIntStand=T, 
                        internalStandID=c("L1","H1","LW"),
                        dataName=dataName,
                        maxPkNum=18, 
                        expectedNonSampPks=7,
                        sdCrefIso.thresh=0.1,
                        sdOrefIso.thresh=0.1,
                        checkRelDiffIntensity=T, # flag for refSamp intensity check
                        refSamp_relDiffInt.thresh=0.75, #*new
                        amplName="Ampl44",
                        ref_relDiffInt.thresh=0.1,
                        sdCsampIso.thresh=0.3, # change for biotic=0.6
                        sdOsampIso.thresh=0.2, # change for biotic=0.6
                        outPath=resultsPath, # *new, specify path to results directory to write files
                        verbose=T) 
end<-Sys.time()
QAtime<-end-start
QAtime # can be seconds to several minutes depending on directory size


# get reference peak data that passed QA/QC
refs.dat<-currFiltered[[1]]
dim(refs.dat)
#[1] 2175   44

# get sample peak data 
samps.dat<-currFiltered[[2]]
dim(samps.dat)
#[1] 3915   44
colnames(samps.dat)


### internal standards check results
# calibration: currently after QA/QC, will include as an option in future version
int_stand.list <- currFiltered[[3]]
# internal standards for each sample peak in an experiment
int_stand.list[[1]][[1]]
#   Analysis Identifier1 PkNr          d18O16O           d13C12C
# 1     2118          H1    7 4.63408927757558 -11.6554800648649
# 2     2118          H1    8 4.57436254703358 -11.6447678139653
# 3     2118          H1    9 4.61171923077752 -11.6692981344158
# 4     2118          H1   10 4.64576054221344 -11.6502162196429
# 5     2118          H1   11 4.63423343292502 -11.6225052367925
# 6     2118          H1   12 4.57530067636935  -11.664273662242
# 7     2118          H1   13  4.6707544365685 -11.5867468660957
# 8     2118          H1   14 4.59079906736126 -11.6385947049331
# 9     2118          H1   15 4.60301367958449 -11.5811913545488
length(int_stand.list[[1]]) # number of experiments that are internal standards
# [1] 79

# average delta for internal standards
int_stand.list[[2]]
#   standard accepted_d18O16O measured_d18O16O
# 1       L1            -8.55        -9.135153
# 2       H1             4.85         4.213694
# 3       LW            -3.85        -4.628521

# delta SD for internal standards
int_stand.list[[3]]
#   standard accepted_SD_d18O16O calculated_SD_d18O16O
# 1       L1                 0.2            0.03996445
# 2       H1                 0.2            0.04171608
# 3       LW                 0.2            0.04295242

# linear model for calibration
standlm<-int_stand.list[[4]]
standlm[[1]]
# Coefficients:
# (Intercept)  measured_d18O16O  
# 0.6701            1.0011  

# stats for linear model
standlm[[2]]
# Residuals:
# 1        2        3 
# -0.07512 -0.03829  0.11341 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)       0.67007    0.09408   7.122  0.08880 . 
# measured_d18O16O  1.00107    0.01472  68.029  0.00936 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1413 on 1 degrees of freedom
# Multiple R-squared:  0.9998,	Adjusted R-squared:  0.9996 
# F-statistic:  4628 on 1 and 1 DF,  p-value: 0.009357



##### visualize pass/fail in pdfs of chromatogram thumbnails
# Get all filenames of .dxf files in the directory
passFiles <- unique(samps.dat$fileId)

# get passed and failed files from previously defined file names
passedInd <- which(fileNames %in% passFiles)
failedFiles <- fileNames[-passedInd]

# get intensity vs time data for passed and failed experiments
passed_raw.list <- raw_data_all(passFiles, path=dataPath)
failed_raw.list <- raw_data_all(failedFiles, path=dataPath)

# plot all raw data and write to pdfs
# passed QC
generic_plot_all_raw(raw.list=passed_raw.list, path=plotsPath,
                     write_pdf=T,
                     pdf_name="passed_qc_spectra.pdf") 
# failed QC
generic_plot_all_raw(failed_raw.list, path=plotsPath,
                     write_pdf=T,
                     pdf_name="failed_qc_spectra.pdf") 



