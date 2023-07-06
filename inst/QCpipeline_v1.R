# dependencies 
#devtools::install_github("isoverse/isoreader")
library(isoreader)
library(dplyr)
library(pracma)

# set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load functions
#source("QAQC_IRMS_functions_v1.R")

##################

### basics - working with one dxf file 

#dataPath <- "./data/dxf_files/abiotic/"
dxf_file <- "170506_NaHCO3 L + NaCl U_.dxf" 
# sorry about the file names, done before me
dataPath <- system.file("extdata/dxf_files/abiotic", dxf_file, package = "QCIRMS")

#file_path <- paste(dataPath,dxf_file,sep="")

# Can print summaries of any of the data files - info on contents
file.summ<-read_summary(dataPath)
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
head(vend.df)[1:3,1:10]
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


# Can plot raw data (Intensity (mV) vs Rt)
generic_raw_plot(raw.df,file_path)
## save to pdf
setwd("./plots/")
pdf(file="one_dxf.pdf",width=6,height=4)
generic_raw_plot(raw.df,file_path)
dev.off()


# Can plot raw data of all files in a vector or directory
# first get all filenames of .dxf files in the directory
# path to dxf files (raw IRMS experimental data)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(dataPath)
fileNames<-all_dxf_files() #uses current working directory
rawList<-raw_data_all(fileNames)

# plot all raw data using code that writes pdfs
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./plots/")
pdf(file="all_dxf.pdf",width=6,height=4)
generic_plot_all_raw(rawList) 
dev.off()



###################
### QA/QC for a directory of volatile CO2 IRMS experiments

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# get expected times for reference peaks
expRef.df<-read.table("./data/qc/referencePeaks_expectedTimes.txt")
expRef.df
#   Ref_Peak_Nr Expected_Start Expected_Rt Expected_End
# 1           1        27.1700     47.3888      50.4929
# 2           2        67.0193     87.1840      90.3422
# 3           4       166.5653    186.7608     189.8572
# 4           5       206.4920    226.5096     229.6987
# 5          16       843.3227    863.5029     866.6069


# provide dataset name and set working directory to data location
dataName<-"abiotic"
setwd(dataPath)

### run QA/QC on directory of dxf files 
# output written to directory where the data is
# execute start <- ... to QAtime together (select then run) to print runtime
start<-Sys.time()
currFiltered<-QAQC_IRMS(unfilteredPath=dataPath,
                        expRef.df=expRef.df, 
                        checkIntStand=F, #internalStandID=c("L1","H1","LW"),
                        dataName=dataName,
                        maxPkNum=18, 
                        expectedNonSampPks=7,
                        sdCrefIso.thresh=0.1,
                        sdOrefIso.thresh=0.1,
                        checkRelDiffIntensity=T,
                        amplName="Ampl44",
                        relDiffInt.thresh=0.1,#for reference peaks, ref:samp hard-coded for now
                        sdCsampIso.thresh=0.3,#change for biotic=0.6
                        sdOsampIso.thresh=0.2,# change for biotic=0.6
                        verbose=T)
end<-Sys.time()
QAtime<-end-start
QAtime # can be seconds to several minutes depending on directory size



## visualize pass/fail in pdfs of chromatogram thumbnails
# Get all filenames of .dxf files in the directory
refs.dat<-currFiltered[[1]]
dim(refs.dat)
#[1] 2175   44

samps.dat<-currFiltered[[2]]
dim(samps.dat)
#[1] 3915   44
colnames(samps.dat)

# get file names of passed and failed samples for plots or further analysis
passFiles <- unique(samps.dat$fileId)

# get failed
dxfFileNames<-all_dxf_files()
passedInd <- which(dxfFileNames %in% passFiles)
failedFiles <- dxfFileNames[-passedInd]

passed_raw.list<-raw_data_all(passFiles)
failed_raw.list<-raw_data_all(failedFiles)

# plot all raw data and write to pdfs
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./plots/")
pdf(file="passed.pdf",width=6,height=4)
generic_plot_all_raw(passed_raw.list) 
dev.off()

pdf(file="failed.pdf",width=6,height=4)
generic_plot_all_raw(passed_raw.list) 
dev.off()

