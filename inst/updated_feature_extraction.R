library(dplyr)
library(isoreader)
library(tsfeatures)
library(doParallel)
library(MASS)
library(ramify)
source("feature_extraction_functions.R")

# set working directory - sets wd of this script wherever you saved it 
setwd("/home/lily/Desktop/EuropaMLMS/SU23/feature_extraction/") #USER EDIT:

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
  #curr_dir
  ## Read in data 
  fileName <- paste(curr_dir, "_samps_allChecks.csv" , sep="")#USER EDIT
  #fileName
  oldwd <- getwd()
  setwd(paste("./QCresults/", curr_dir , sep=""))
  newdata<-read.csv(fileName)
  # get Analysis numbers
  analysisNums.vec<-unique(newdata$Analysis)
  #length(analysisNums.vec)
  #dim(newdata)
  setwd(oldwd)
  #getwd()
  # need to find the original dxf files for each analysis for raw data
  dxfDir<-paste(getwd(),"/new_ocean_worlds/", curr_dir ,sep="") #USER EDIT
  #dxfDir
  #setwd(dxfDir)
  
  ##### Pre-extraction processing
  # separate data by experiment
  fullData.list<-separate_by_analysis(newdata)
  #length(fullData.list)
  # get raw data for feature extraction
  rawTsDat.list<-rawDat_by_analysis_num(analysisNums.vec,dxfDir)
  #length(rawTsDat.list)
  
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
  #extrFeat.df
  #getwd()
  # add to output list
  feat.list[[j]]<-extrFeat.df
  # write to file
  writeName <- paste("./feature_extraction_results/", curr_dir,"_extracted_ts_feat.csv", sep="") #USER EDIT
  #writeName
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
