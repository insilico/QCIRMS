library(dplyr)
library(isoreader)
library(tsfeatures)
library(doParallel)
library(MASS)
library(ramify)
source("feature_extraction_functions.R")

# set working directory - sets wd of this script wherever you saved it 
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/lily/Desktop/EuropaMLMS/SU23/feature_extraction/")
## vector of directories that could be calibrated 
dxf_dirs.df <- read.csv("./QCresults/new_OW_passed_dirs.csv")
dxf_dirs.vec <- dxf_dirs.df$passed_dirs
cat("numbers of directories for feature extraction: ",length(dxf_dirs.vec))

# set up the parallelization
numCores <- detectCores()-2
cat("Number of available cores:", numCores) #40
#print(paste("Using half the number of available cores: ",numCores,sep=""))
registerDoParallel(numCores)

print("Beginning feature extraction...")
startTime<-Sys.time()

for(i in seq(1, length(dxf_dirs.vec))){
  curr_dir <- dxf_dirs.vec[i] 
  curr_dir
  ## Read in data 
  fileName <- paste(curr_dir, "_samps_allChecks.csv" , sep="")
  fileName
  oldwd <- getwd()
  setwd(paste("./QCresults/", curr_dir , sep=""))
  newdata<-read.csv(fileName)
  # get Analysis numbers
  analysisNums.vec<-unique(newdata$Analysis)
  #length(analysisNums.vec)
  dim(newdata)
  setwd(oldwd)
  getwd()
  # need to find the original dxf files for each analysis for raw data
  dxfDir<-paste(getwd(),"/new_ocean_worlds/", curr_dir ,sep="") #LAC - needs to be edited
  dxfDir
  #setwd(dxfDir)
  
  ##### Pre-extraction processing
  # separate data by experiment
  fullData.list<-separate_by_analysis(newdata)
  length(fullData.list)
  # get raw data for feature extraction
  rawTsDat.list<-rawDat_by_analysis_num(analysisNums.vec,dxfDir)
  length(rawTsDat.list)
  
  ##### Feature extraction
  feat.list<-list()
  # initialize the ts output file
  j=1
  intDat44<-rawTsDat.list[[j]]
  #head(intDat44)
  analysisNum<-colnames(rawTsDat.list[[j]])
  # extract features
  extrFeat.df<-extract_tsfeatures(intDat44)
  extrFeat.df$Analysis<-analysisNum
  extrFeat.df
  #getwd()
  # add to output list
  feat.list[[j]]<-extrFeat.df
  # write to file
  writeName <- paste("./feature_extraction_results/", curr_dir,"_extracted_ts_feat.csv", sep="")
  writeName
  write.table(extrFeat.df, writeName,
              col.names=T, quote=F, row.names=F, sep=",")
  print(paste("Successful extraction of ",j,"/",length(rawTsDat.list),
              " experiments in directory ", i, "/", length(dxf_dirs.vec), 
              "... ", sep=""))
  # now run parallel feature extraction for the rest of the data
  foreach(j=seq(2,length(rawTsDat.list)))%dopar%{ 
    #for(i in seq(2,5)){
    # get intensity data
    intDat44<-rawTsDat.list[[j]]
    analysisNum<-colnames(rawTsDat.list[[j]])
    # extract features
    extrFeat.df<-extract_tsfeatures(intDat44)
    extrFeat.df$Analysis<-analysisNum
    # add to output list
    feat.list[[j]]<-extrFeat.df
    # write
    write.table(extrFeat.df, writeName, append=T, col.names=F, quote=F, row.names=F, sep=",")
    print(paste("Successful extraction of ", j, "/", length(rawTsDat.list),
                " experiments in directory ", i, "/", length(dxf_dirs.vec), 
                "... ", sep=""))
  }
  
}
endTime<-Sys.time()
timeDiff<-endTime-startTime
print(paste("feature extraction: ",timeDiff,sep=""))

##########
# MS features - avg and sd of specified features
# this code uses calibrated isotopes -- may not be available (edit fileName)

for(i in seq(1, length(dxf_dirs.vec))){
  curr_dir <- dxf_dirs.vec[i]
  ## Read in data
  fileName <- paste(curr_dir, "_calibr_samps_allChecks.csv" , sep="")
  oldwd <- getwd()
  setwd(paste("./QCresults/", curr_dir , sep=""))
  newdata<-read.csv(fileName)
  # get Analysis numbers
  analysisNums.vec<-unique(newdata$Analysis)
  setwd(oldwd)

  # separate data by experiment
  fullData.list<-separate_by_analysis(newdata)
  length(fullData.list)
  # colnames to not take avg and sd of
  rm_cols <- c("fileId", "Identifier1", "Analysis", "Preparation",
               "DateTime", "PeakNr", "Start", "Rt", "End", "ListFirstPeak",
               "IsRef", "RefName", "Rps45CO244CO2", "Rps46CO244CO2")
  add_back <- c("fileId", "Identifier1", "Analysis", "Preparation", "DateTime")
  red_meas.list <- list()
  # loop though the dataframes in list
  for(j in seq(1,length(fullData.list))){
    full_dat <- fullData.list[[j]]
    # take averages and sd of cols
    curr_analysis <- curr_dat$Analysis[1]
    curr_dat <- full_dat[,-which(colnames(full_dat) %in% rm_cols)]
    new_df <- as.data.frame(matrix(rep(NA,2*dim(curr_dat)[2]),ncol=2*dim(curr_dat)[2]))
    colnames.vec <- c()
    vals.vec <- c()
    # loop through cols
    for(k in seq(1,length(colnames(curr_dat)))){
      curr_col <- colnames(curr_dat)[k]
      curr_col.vec <- curr_dat[,which(colnames(curr_dat)==curr_col)]  
      avg_col <- mean(curr_col.vec)
      vals.vec <- c(vals.vec, avg_col)
      avg_name <- paste("avg_", curr_col, sep="")
      colnames.vec <- c(colnames.vec,avg_name)
      sd_col <- sd(curr_col.vec)
      vals.vec <- c(vals.vec, sd_col)
      sd_name <- paste("sd_", curr_col, sep="")
      colnames.vec <- c(colnames.vec,sd_name)
    }
    # add colnames + data, add df to list
    add_back.dat <- full_dat[1,which(colnames(full_dat) %in% add_back)]
    add_back.dat
    new_df[1,] <- vals.vec
    colnames(new_df) <- colnames.vec
    new_full.df <- cbind.data.frame(add_back.dat,new_df)
    new_full.df
    red_meas.list[[j]] <- new_full.df
  }
  # list is experiments in dxf dir
  # list to df
  red_meas.df <- do.call("rbind",red_meas.list)
  # write file in feature_extraction directory
  outFileName <- paste("./feature_extraction/", curr_dir, "_ms_feat.csv", sep="")
  outFileName
  write.table(red_meas.df, outFileName, row.names=F, quote=F, sep=",")
}
##############