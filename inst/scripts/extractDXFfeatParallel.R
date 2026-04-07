

library(tsfeatures)
library(foreach)
library(doParallel)
library(MASS)
library(ramify)
library(isoreader)

# set up the parallelization
numCores <- detectCores()#/2
print(paste("Number of cores detected:",numCores),sep="") #40
#print(paste("Using half the number of available cores: ",numCores,sep=""))
registerDoParallel(numCores)

## **EDIT:filenames and varibales
## Read in data 
newdata<-read.csv("bioAbio_samps_QC_2022-06-25.csv",header=T)
#head(newdata)
# need to find the original dxf files for each analysis for raw data
dxfDir<-"~/Desktop/EuropaMLMS/microbe/QC_SU22/bio_abio/bioAbio/bioAbioDXF"

# variables to use for feature extraction
#vars.vec<-c("d13C12C","d18O16O")
#vars.vec<-c("d13C12C","d18O16O","d18O13C")
#rmExtraVars.vec<-c("PeakNr") #"d18O13C" # bc making the data 1-D

# get analysis numbers to grab raw intensity vs time data
analysisNums.vec<-unique(newdata$Analysis)
#head(analysisNums.vec)

# read in avg delta values to build final ML dataframe
#avgDeltas.df<-read.csv("avgDeltasDXF.csv",header=T)
#head(avgDeltas.df)
# Analysis avg_d18O16O avg_d13C12C avg_d18O13C
#1     2739   -7.644012   -11.65822   0.6556780
#2     2742   -7.678239   -11.66683   0.6581270
#3     2758   -7.295827   -11.53062   0.6327352
#4     3240   -6.643230   -11.92245   0.5572045
#5     3250   -6.582490   -11.97722   0.5495850
#6     3263   -6.565071   -11.91321   0.5510763
#avgDeltas.df<-read.csv("microbeAvgDeltasData.csv",header=T)
#head(avgDeltas.df)
# Analysis avg_d18O16O avg_d13C12C avg_d18O13C
#1     2739   -7.644012   -11.65822   0.6556780
#2     2742   -7.678239   -11.66683   0.6581270
#3     2758   -7.295827   -11.53062   0.6327352
#4     3240   -6.643230   -11.92245   0.5572045
#5     3250   -6.582490   -11.97722   0.5495850
#6     3263   -6.565071   -11.91321   0.5510763
#dim(avgDeltas.df)
#avgAmpls.df<-read.csv("microbeAvgAmplData.csv",header=T)
#head(avgAmpls.df)
# Analysis avg_Ampl44 avg_Ampl45 avg_Ampl46
#1     2739   3659.850   4312.484   5153.190
#2     2742   3657.459   4309.485   5149.380
#3     2758   3696.565   4353.240   5204.446
#4     3240   2518.434   2967.673   3551.860
#5     3250   3714.100   4376.742   5238.228
#6     3263   3694.329   4356.123   5214.847
#dim(avgAmpls.df)

##### Functions --
# need to make a list for tsfeatures; separate data by experiment
separate_by_analysis<-function(raw.df){
  # make a list of experimental data separated by analysis numbers
  analyses<-raw.df$Analysis
  totNumRows<-length(analyses)
  differentAnalyses<-unique(analyses)
  analysisData.list<-list()
  analysisSetInd.vec<-c()
  uniqueAnalysisInd<-1
  # check for end of analysis numbers
  for(i in seq(2,(totNumRows+1))){
    prevAnalysis<-analyses[i-1]
    currAnalysis<-analyses[i]
    # check for end of analysis nums
    if(i==(totNumRows+1)){
      #add last element to vec and then to list
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
      analysisData.list[[uniqueAnalysisInd]]<-raw.df[analysisSetInd.vec,]
      break
    }
    # check if the analysis numbers are the same
    if(currAnalysis==prevAnalysis){
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
    } else{ # different analysis numbers
      # add i-1 to vec
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
      # add analysis set to list
      analysisData.list[[uniqueAnalysisInd]]<-raw.df[analysisSetInd.vec,]
      # reset analysis index vector and set next list index
      analysisSetInd.vec<-c()
      uniqueAnalysisInd<-uniqueAnalysisInd+1
    }
  }
  return(analysisData.list)
}
# get iso file info for vec of files
file_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file_info<-msdat %>%
    iso_get_file_info(
      select = c(
        #rename?
        Identifier_1 = `Identifier 1`,
        # select columns without renaming
        `Analysis`, `Preparation`,
        # select the time stamp and rename it to `Date & Time`
        Date_and_Time = file_datetime
      ),
      # explicitly allow for file specific rename (for the new ID column)
      file_specific = TRUE #mostly useful with data from different instruments
    )
  # convert from tibble to df
  file_info.df<-as.data.frame(file_info)
}
# get raw intensity vs time for vec of files
raw_data_all<-function(files){
  raw.list<-list()
  for(i in seq(1,length(files))){
    raw.dat<-raw_data(files[i])
    raw.list[[i]]<-raw.dat
  }
  return(raw.list)
}
# for one file
raw_data<-function(file){
  #num_files<-length(files)
  msdat<-iso_read_continuous_flow(file)#files[1:num_files]
  raw_dat<- msdat %>% iso_get_raw_data()
  raw_dat.df<-as.data.frame(raw_dat)
}
# get the draw data for vec of analyses
rawDat_by_analysis_num<-function(analysisNums.vec,dxfDir){
  # get files in the dxfDir
  fileNames<-list.files(dxfDir)
  dxfFilesInd<-which(grepl(".dxf",fileNames))
  dxfFiles<-fileNames[dxfFilesInd]
  numFiles<-length(dxfFiles)
  # wd
  oldDir<-getwd()
  setwd(dxfDir)
  #
  numAnalyses<-length(analysisNums.vec)
  # use isoreader to get file info so can search for analysis
  dxfFilesInfo.df<-file_info(dxfFiles)
  #head(dxfFilesInfo.df)
  # get all analysis numbers
  dxfAnalysis.vec<-dxfFilesInfo.df$Analysis
  dxfFiles<-dxfFilesInfo.df$file_id
  # search for analysis numbers match -- hopefully find 195
  analysisFileNames.vec<-c()
  matchAnalysis.vec<-c()
  for(i in seq(1,numAnalyses)){
    if(analysisNums.vec[i] %in% dxfAnalysis.vec){
      matchInd<-match(analysisNums.vec[i],dxfAnalysis.vec)
      matchAnalysis.vec<-c(matchAnalysis.vec,analysisNums.vec[i])
      # get filename
      fileMatch<-dxfFiles[matchInd]
      analysisFileNames.vec<-c(analysisFileNames.vec,fileMatch)
    }
  }
  
  # get raw ms data for analysis we have
  analysisRawData.list<-raw_data_all(analysisFileNames.vec)  
  
  # get intensity data 
  intDat.list<-list()
  fileID.vec<-c()
  for(i in seq(1,length(analysisRawData.list))){
    intDat44<-analysisRawData.list[[i]]$v44.mV
    intDat44.df<-as.data.frame(matrix(intDat44,ncol=1))
    # get the file name too
    file_id<-analysisRawData.list[[i]]$file_id[1]
    fileID.vec<-c(fileID.vec,file_id)
    # add analysis number!
    anNum<-matchAnalysis.vec[i]
    colnames(intDat44.df)<-anNum
    # add to list
    intDat.list[[i]]<-intDat44.df
  }
  # reset working directory
  setwd(oldDir)
  return(intDat.list)
}
# Feature extraction
extract_tsfeatures<-function(intensity.data){
  features.tib<-tsfeatures(intensity.data,
                           features=c("acf_features","arch_stat","crossing_points",
                                      "entropy","flat_spots","heterogeneity",
                                      "holt_parameters","hurst",
                                      "lumpiness","max_kl_shift","max_level_shift",
                                      "max_var_shift","nonlinearity","pacf_features",
                                      "stability",
                                      "stl_features",
                                      "unitroot_kpss",
                                      "unitroot_pp",
                                      "ac_9",
                                      "firstmin_ac",
                                      "firstzero_ac",
                                      "fluctanal_prop_r1",
                                      "histogram_mode","localsimple_taures","motiftwo_entro3",
                                      "outlierinclude_mdrmd","sampenc","sampen_first",
                                      "std1st_der","trev_num","spreadrandomlocal_meantaul",
                                      "walker_propcross"))
  features.df<-as.data.frame(features.tib)
}
##### 

##### Pre-extraction processing
# separate data by experiment
fullData.list<-separate_by_analysis(newdata)

# get raw data for feature extraction
rawTsDat.list<-rawDat_by_analysis_num(analysisNums.vec,dxfDir)

##### Feature extraction
print("Beginning feature extraction...")
startTime<-Sys.time()
feat.list<-list()
# initialize the ts output file
i=1
intDat44<-rawTsDat.list[[i]]
#head(intDat44)
analysisNum<-colnames(rawTsDat.list[[i]])
# extract features
extrFeat.df<-extract_tsfeatures(intDat44)
extrFeat.df$Analysis<-analysisNum
#extrFeat.df
# add to output list
feat.list[[i]]<-extrFeat.df
# write to file
write.table(extrFeat.df,"dxfTStestPar.csv",col.names=T,quote=F,row.names=F,sep=",")
print(paste("Successful extraction of ",i,"/",length(rawTsDat.list),
            " experiments...",sep=""))
# now run parallel feature extraction for the rest of the data
foreach(i=seq(2,length(rawTsDat.list)))%dopar%{ 
  #for(i in seq(2,5)){
  # get intensity data
  intDat44<-rawTsDat.list[[i]]
  analysisNum<-colnames(rawTsDat.list[[i]])
  # extract features
  extrFeat.df<-extract_tsfeatures(intDat44)
  extrFeat.df$Analysis<-analysisNum
  # add to output list
  feat.list[[i]]<-extrFeat.df
  # write
  write.table(extrFeat.df,"dxfTStestPar.csv",append=T,col.names=F,quote=F,row.names=F,sep=",")
  print(paste("Successful extraction of ",i,"/",length(rawTsDat.list),
              " experiments...",sep=""))
}
endTime<-Sys.time()
timeDiff<-endTime-startTime
print(paste("feature extraction: ",timeDiff,sep=""))

# match the extracted data with the original ML labels
# read in table
#extrFeat.tab<-read.csv("dxfTStestPar.csv",header=T)
#head(extrFeat.tab)
# get analysis numbers
#fullDataAn.vec<-c()
#for(i in seq(1,length(fullData.list))){
#  fullDataAn.vec<-c(fullDataAn.vec,fullData.list[[i]]$Analysis[1])
#}
#extrFeatAn.vec<-extrFeat.tab$Analysis
## match analysis numbers and create df for each analysis
# add to final dataframe
#colnames(newdata)

#colnames(extrFeat.tab) #includes Analysis
#final.df<-extrFeat.tab # build on tsfeatures df

#colsToAdd<-c("biotic","fileId","Identifier1",    
#             "Preparation","MgSO4","NaCl",         
##             "NaHCO3","Na2SO4","MgSO4_L","NaCl_U",       
#             "NaCl_L","NaHCO3_L","NaHCO3_U","Na2SO4_L",     
#             "Na2SO4_U","rel_MgSO4","rel_NaCl","rel_NaHCO3",   
#             "rel_Na2SO4","saltContent","anion","cation")
#colsToAdd<-c("fileId","Identifier1","DateTime",    
#  "Preparation","biotic","avg_d13C12C",
#  "avg_d18O13C","avg_d18O16O","avg_Ampl46","avg_IntensityAll",
#  "avg_rIntensityAll","avg_pkArea")


#head(avgDeltas.df)
#deltaColsToAdd<-c("avg_d18O16O", "avg_d13C12C","avg_d18O13C")
#deltaAn.vec<-avgDeltas.df$Analysis
#amplAn.vec<-avgAmpls.df$Analysis
#combDF.list<-list()
#for(i in seq(1,length(extrFeatAn.vec))){
#  
#  matchInd<-which(fullDataAn.vec==extrFeatAn.vec[i])
#  
#  deltaMatchInd<-which(deltaAn.vec==extrFeatAn.vec[i])
#  amplMatchInd<-which(amplAn.vec==extrFeatAn.vec[i])
#  
#  match.df<-fullData.list[[matchInd]]
#  keepInd<-which(colnames(match.df) %in% colsToAdd)
#  match.df<-match.df[1,keepInd]
#  #match.df
  
  # rm analysis from deltas and ampl bc its in tsfeatures
#  deltaAnInd<-which(colnames(avgDeltas.df)=="Analysis")
#  deltaMatch.df<-avgDeltas.df[deltaMatchInd,]
#  deltaMatch.df<-deltaMatch.df[,-deltaAnInd]
  #deltaMatch.df
  ##
#  amplAnInd<-which(colnames(avgAmpls.df)=="Analysis")
#  amplMatch.df<-avgAmpls.df[amplMatchInd,]
#  amplMatch.df<-amplMatch.df[,-amplAnInd]
#  #amplMatch.df
  # construct final df
#  final.df<-as.data.frame(c(extrFeat.tab[i,],
#                            deltaMatch.df, amplMatch.df,
#                            match.df))

  #final.df
 # combDF.list[[i]]<-final.df
#}

# flatten list to df
#for(i in seq(1,length(combDF.list))){
##  if(i==1){
#    write.table(combDF.list[[i]],"dxfTSsemiFinal.csv",col.names=T,quote=F,row.names=F,sep=",")
#  }else{
#    write.table(combDF.list[[i]],"dxfTSsemiFinal.csv",append=T,col.names=F,quote=F,row.names=F,sep=",")
#  }
#}
#semiFinal.df<-read.csv("dxfTSsemiFinal.csv",header=T)
#colnames(semiFinal.df)
#origColNames<-colnames(semiFinal.df)

# ** 
#semiFinal.df<-final.df
#origColNames<-colnames(final.df)
# **

# organize df
#colOrder<-c(# labels:
#            "fileId", "DateTime","Analysis",
#            "Identifier1",
#            "Preparation",
#            "biotic",
#            
            #"NaHCO3","NaCl","Na2SO4","MgSO4",#"KCl","MgCl2",
            #"NaHCO3_L","NaHCO3_U","NaCl_L","NaCl_U",
            #"Na2SO4_L","Na2SO4_U","MgSO4_L",#"KCl_L","KCl_U",
            #"MgCl2_L","MgCl2_U",
            #"rel_NaHCO3","rel_NaCl","rel_Na2SO4","rel_MgSO4",
            #"rel_KCl","rel_MgCl2", 
            #"saltContent","anion","cation",
            
            # numerical data
            # new calculated data from sample peaks
 #           "avg_d18O16O","avg_d13C12C","avg_d18O13C",
            #"avg_Ampl44","avg_Ampl45","avg_Ampl46",
 #           "avg_d13C12C",
#            "avg_d18O13C","avg_d18O16O","avg_Ampl46","avg_IntensityAll",
#            "avg_rIntensityAll","avg_pkArea",
            
            # time series nuemrical data
#            "x_acf1","x_acf10","diff1_acf1","diff1_acf10",
#            "diff2_acf1","diff2_acf10","ARCH.LM","crossing_points",
##            "entropy","flat_spots","arch_acf","garch_acf","arch_r2",
#            "garch_r2","alpha","beta","hurst","lumpiness",
#            "max_kl_shift", "time_kl_shift", "max_level_shift",
#            "time_level_shift", "max_var_shift","time_var_shift","nonlinearity",
#            "x_pacf5","diff1x_pacf5","diff2x_pacf5","stability","nperiods",
#            "seasonal_period","trend","spike","linearity","curvature",
#            "e_acf1","e_acf10","unitroot_kpss","unitroot_pp","ac_9",
 #           "firstmin_ac","firstzero_ac","fluctanal_prop_r1","histogram_mode",
#            "localsimple_taures","motiftwo_entro3","outlierinclude_mdrmd",
##            "sampenc", "sampen_first", "std1st_der","trev_num",
#            "spreadrandomlocal_meantaul","walker_propcross"
#)
#finalOrg.df<-as.data.frame(matrix(rep(NA,dim(semiFinal.df)[1]*dim(semiFinal.df)[2]),
#                                  ncol=dim(semiFinal.df)[2]))
#colnames(finalOrg.df)<-colOrder
#for(i in seq(1,length(colOrder))){
#  orderInd<-which(origColNames==colOrder[i])
#  currCol<-semiFinal.df[,orderInd]
#  finalOrg.df[,i]<-currCol
#}
#finalOrg.df

# write final dataframe to file
#write.table(finalOrg.df,"dxfTSfeatFullFinal.csv",row.names=F,quote=F,sep=",")
###
