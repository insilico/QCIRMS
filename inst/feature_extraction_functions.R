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