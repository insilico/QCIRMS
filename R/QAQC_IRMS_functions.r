############
# (1)
#' combineVendFileInfo: function that processes vendor data for a directory of dxf files
#' @param path string that is the path to the directory of .dxf files
#' @param combColNames vector of combined column names for which to extract data from the dxf vendor table
#' @param outputID string that identifies the dataset 
#' @return list of IRMS analysis vendor table dataframes 
#' @examples 
#' Usage example
#' vend.list <- combineVendFileInfo(path = "./abiotic/", 
#'                                  combColNames = c("fileId","Identifier1","Analysis","Preparation","DateTime",
#'                                   "PeakNr","Start","Rt","End","Ampl44","Ampl45",
#'                                   "Ampl46","BGD44","BGD45","BGD46","rIntensity44","rIntensity45",
#'                                   "rIntensity46","rIntensityAll","Intensity44","Intensity45",
#'                                   "Intensity46","IntensityAll","ListFirstPeak","rR45CO244CO2",
#'                                   "rR46CO244CO2","IsRef","R45CO244CO2","RefName","rd45CO244CO2",
#'                                   "d45CO244CO2", "R46CO244CO2", "rd46CO244CO2","d46CO244CO2",
#'                                   "R13C12C","d13C12C","AT13C12C","R18O16O","d18O16O", "AT18O16O",
#'                                   "R17O16O","d17O16O","Rps45CO244CO2","Rps46CO244CO2"),
#'                                  outputID = "abiotic")
#' @export
combineVendFileInfo<-function(path,combColNames,outputID){
  # make a combined dataframe for vendor and file info
  files<-all_dxf_files()
  vend<-vendor_info_all(files)
  fileInfo<-file_info(files)
  # make a list of combined vendor data and file info dfs for each experiment
  combList<-list()
  failedInd<-c()
  failedFiles<-c()
  for(i in seq(1,length(files))){
    # make file info df to add to vend df
    mat<-matrix(rep(NA,nrow(vend[[i]])*ncol(fileInfo)),#ncol(fileInfo[[i]])),
                nrow=nrow(vend[[i]]))
    df<-as.data.frame(mat)
    df[1:nrow(df),]<-fileInfo[i,]
    # need to check if there is data
    isVendData<-!(dim(vend[[i]])[1]==0)
    if(isVendData){
      # make a combined df
      combDF<-data.frame(cbind(df,vend[[i]][,2:ncol(vend[[i]])]))
      # fileId, Identifier1, Analysis, DateTime, PeakNr,
      # Start, Rt, End, Ampl44,
      colnames(combDF)<-combColNames
      # add to list
      combList[[i]]<-combDF
    }else{
      fileName<-files[i]
      print(paste("No data found in file: ",fileName,". File name inserted at index ",i,sep=""))
      failedFiles<-c(failedFiles,fileName)
      failedInd<-c(failedInd,i)
      combList[[i]]<-fileName
    }
  }
  # if there is missing data remove that file
  if(length(failedInd)>0){
    print("Removing files with no data from the analysis...")
    onlyDat<-combList[c(-failedInd)]
    failedFiles.df<-as.data.frame(matrix(failedFiles),ncol=1)
    colnames(failedFiles.df)<-c("No_Data")
    write.table(failedFiles.df,file=paste(outputID,"_NoData.csv",sep=""),row.names=F,quote=F,sep=",")
    return(onlyDat)
  }else{
    return(combList)
  }
}


# (2)
#' removeFailedAnalysesDXF: function that removes experiments that fail QA/QC from the processing pipeline
#' @param sepList list of IRMS vendor dataframes separated by experiment
#' @param expRef.df dataframe of expected reference peak times with column names = c("Ref_Peak_Nr", "Expected_Start", "Expected_Rt", "Expected_End")
#' @param maxPkNum integer for the maximum number of acceptable peaks in a spectrum; default value = 18
#' @param expRefPkNr integer for the expected number of reference peaks; default value = 5
#' @param diff.t the accepted interval of peak retention times in seconds; default vale = 10
#' @param sdCrefIso.thresh the accepted value for reference peak delta 13C ; default value = 0.1
#' @param expectedNonSampPks the expected number of non-sample peaks; default value = 7
#' @param sdOrefIso.thresh the accepted value for reference peak delta 18O; default value = 0.1
#' @param relDiffInt.thresh the accepted relative difference in peak intensity; default value = 0.1 (10% difference max tolerated)
#' @param amplName string for the column to use for the intensity check; default value = "Ampl44"
#' @param sdCsampIso.thresh the accepted value for sample peak delta 13C; default value = 0.3
#' @param sdOsampIso.thresh the accepted value for sample peak delta 18O; default value = 0.2
#' @param flushExpT the expected retention time in seconds for the flush decontamination peak, if present; default value = 135
#' @param flushTint the accepted interval in seconds for the flush peak retention time; default value = 15
#' @param firstSampExpT the expected retention time for the first sample peak in seconds; default value = 275
#' @param firstSampTint the accepted interval in seconds for the arrival of the first sample peak; default value = 15
#' @param verbose whether to print information about failed analyses; default value = T
#' @return list of eight elements containing information on samples that failed the checks and the data that passed
#' ret.list[[1]]: info for samples that failed the peak number check
#' ret.list[[2]]: info for samples that failed the reference peaks check
#' ret.list[[3]]: info for samples that failed the reference isotope standard deviation check
#' ret.list[[4]]: info for samples that failed the reference peak relative intensity check 
#' ret.list[[5]]: info for samples that failed the sample peak processing steps 
#' ret.list[[6]]: info for samples that failed the sample peak isotope standard deviation check
#' ret.list[[7]]: list of reference peak data that passed the six checks above
#' ret.list[[8]]: list of sample peak data that passed the six checks above
#' @examples 
#' Usage example
#' processed.list <- removeFailedAnalysesDXF(sepList=combList,
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
#'                                           verbose=T)
#' dont export yet
removeFailedAnalysesDXF<-function(sepList, 
                                  expRef.df, 
                                  maxPkNum=18, 
                                  expRefPkNr=5,
                                  diff.t=10,
                                  sdCrefIso.thresh=0.1, 
                                  expectedNonSampPks=7,
                                  sdOrefIso.thresh=0.1,
                                  relDiffInt.thresh=0.1,
                                  amplName="Ampl44",
                                  sdCsampIso.thresh=0.3,
                                  sdOsampIso.thresh=0.2,
                                  flushExpT=135,
                                  flushTint=15,
                                  firstSampExpT=275,
                                  firstSampTint=15,
                                  verbose=T){
  sep.list<-list()
  samps.list<-list()
  refs.list<-list()
  # get reference expected times
  expRefStart<-expRef.df$Expected_Start
  expRefRt<-expRef.df$Expected_Rt
  expRefEnd<-expRef.df$Expected_End
  # set up loop
  numCheck<-length(sepList)
  # initialize list indices
  refIndex<-1
  sampIndex<-1
  peakIndex<-1
  # intialize vecs for each check-analyis nums and identifier1
  # peak number check-0 peaks and over max
  failedPkNum.vec<-c()
  failedPkID.vec<-c()
  failedNumPks.vec<-c()
  # reference times check
  failedRefAnNum.vec<-c()
  failedRefID.vec<-c()
  failedRefDat.list<-list()
  # reference isotope ratios
  failedRefIso.vec<-c()
  failedRefIsoID.vec<-c()
  failedRefIsoSD.vec<-c() #d18O/16O
  failedRefIsoSDC.vec<-c()
  failedRefIsoReason.vec<-c()
  # reference intensity similiarity
  failedRefInt.vec<-c()
  failedRefIntID.vec<-c()
  failedRefIntRD.vec<-c()
  # sample processing
  failedSampAnNum.vec<-c()
  failedSampID.vec<-c()
  failedSampReas.vec<-c()
  # sample isotope ratio check
  failedSampIso.vec<-c()
  failedSampIsoID.vec<-c()
  failedSampIsoSD.vec<-c()
  failedSampIsoSDC.vec<-c()
  failedSampIsoReason.vec<-c()
  
  # TODO: sampRef peak int check?
  
  for(i in seq(1,numCheck)){
    # check number of peaks
    numPeaksCheck<-maxNumPeaksDXF(sepList[[i]],
                                  maxExpectedPks=maxPkNum,
                                  expectedNonSampPks=expectedNonSampPks,
                                  verbose=verbose)
    numPks<-numPeaksCheck[[2]]$NumPks
    numPksPass<-numPeaksCheck[[1]]
    if(!numPksPass){ 
      failedPkNum.vec<-c(failedPkNum.vec,sepList[[i]]$Analysis[1])
      failedPkID.vec<-c(failedPkID.vec,sepList[[i]]$Identifier1[1])
      failedNumPks.vec<-c(failedNumPks.vec,numPks)
    } else{ # add to data list for further checks
      sep.list[[peakIndex]]<-sepList[[i]]
      peakIndex<-peakIndex+1
    }
  }
  
  # update numCheck after removing bad data
  numCheck<-length(sep.list)
  # check for any data, then proceed with the rest of the checks
  if(length(sep.list)>0){
    if(is.null(failedPkNum.vec)){
      if(verbose==T){print("No analyses failed max total peak number check.")}
      failedPkNum.vec<-c(NA)
      failedPkID.vec<-c(NA)
      failedNumPks.vec<-c(NA)
    }
    # remove spaces in Identifier1 and replace with underscores
    if(sum(grepl(" ", failedPkID.vec))>0){
      failedPkID.vec<-gsub(" ","_",failedPkID.vec)
    }
    failedPkNum.df<-as.data.frame(matrix(c(failedPkNum.vec,failedPkID.vec,failedNumPks.vec),ncol=3))
    colnames(failedPkNum.df)<-c("failed_PkNr_Analysis","failed_PkNr_Identifier1","failed_PkNr_NumPeaks")
    #failedPkNum.df
    
    ## reference and sample peak checks
    failedRefDatIndex<-1
    for(i in seq(1,numCheck)){
      # check reference peaks times
      sampRefCheck<-reference_times_checkDXF(vend.df=sep.list[[i]],
                                             expectedPeak.num=expRefPkNr,
                                             diff.t=diff.t,
                                             expectedStart=expRefStart,
                                             expectedRt=expRefRt,
                                             expectedEnd=expRefEnd,
                                             verbose=verbose)
      if(!sampRefCheck[[1]]){ # if sample failed ref time check
        failedRefAnNum.vec<-c(failedRefAnNum.vec,sep.list[[i]]$Analysis[1])
        failedRefID.vec<-c(failedRefID.vec,sep.list[[i]]$Identifier1[1])
        # remove spaces in Identifier1, replace with underscores
        if(sum(grepl(" ", failedRefID.vec))>0){
          failedRefID.vec<-gsub(" ","_",failedRefID.vec)
        }
        # return failed reference peak data
        failedRefDat.list[[failedRefDatIndex]]<-sampRefCheck[[2]]
        failedRefDatIndex<-failedRefDatIndex+1
      } else{ # perform the rest of the checks
        # add ref iso ratio check and intensity similarity check
        # need pknr.vec for iso_ratio_similarity
        refPkNr.vec<-sampRefCheck[[2]]$PkNr
        # reference peak isotope similarity check
        
        refIsoCheck<-isoR_similarityDXF(vend.df=sep.list[[i]],
                                        peakNr.vec=refPkNr.vec,
                                        sdC.thresh=sdCrefIso.thresh,
                                        sdO.thresh=sdOrefIso.thresh,
                                        verbose=verbose)
        # check if failed one of the checks
        if((!refIsoCheck[[1]]$d18O16O)|(!refIsoCheck[[1]]$d13C12C)){
          refReasInd.vec<-which(refIsoCheck[[1]][1,]==FALSE)
          refReason<-colnames(refIsoCheck[[1]])[refReasInd.vec]
          if(length(refReason)==2){
            # both isotope SDs failed
            refReason<-c("d13C12C_d18O16O")
          }
          failedRefIso.vec<-c(failedRefIso.vec,sep.list[[i]]$Analysis[1])
          failedRefIsoID.vec<-c(failedRefIsoID.vec,sep.list[[i]]$Identifier1[1])
          #*
          if(sum(grepl(" ", failedRefIsoID.vec)>0)){
            failedRefIsoID.vec<-gsub(" ","_",failedRefIsoID.vec)
          }
          #*
          failedRefIsoSD.vec<-c(failedRefIsoSD.vec,refIsoCheck[[2]]$SD_d18O16O)
          failedRefIsoSDC.vec<-c(failedRefIsoSDC.vec,refIsoCheck[[2]]$SD_d13C12C)
          failedRefIsoReason.vec<-c(failedRefIsoReason.vec,refReason)
        } else{ # if it passes check intensity similarities
          # check the reference peak isotope similarities
          refIntCheck<-intensity_similarityDXF(vendAmpl=sep.list[[i]]$Ampl44[refPkNr.vec],
                                               amplName=amplName,
                                               peakNr.vec=refPkNr.vec,
                                               relDiffInt.thresh=relDiffInt.thresh,
                                               verbose=verbose)
          
          # if it failed add to output
          if(!refIntCheck[[1]]){
            failedRefInt.vec<-c(failedRefInt.vec,sep.list[[i]]$Analysis[1])
            failedRefIntID.vec<-c(failedRefIntID.vec,sep.list[[i]]$Identifier.1[1])
            #*
            if(sum(grepl(" ", failedRefIntID.vec)>0)){
              failedRefIntID.vec<-gsub(" ","_",failedRefIntID.vec)
            }
            #*
            failedRefIntRD.vec<-c(failedRefIntRD.vec,refIntCheck[[3]])
          }else{# only if passes all above 3 ref checks add to ref data
            refs.list[[refIndex]]<-sampRefCheck[[3]]
            refIndex<-refIndex+1
            # ** only check samples if ref data was added
            # process the sample peaks: remove flush peak and 1st sample peak (contamination)
            
            sampProcess<-sample_peaks_processDXF(refTimesOutput=sampRefCheck,
                                                 vend.df=sep.list[[i]],
                                                 flushExpT=flushExpT, flushTint=flushTint,
                                                 firstSampExpT=firstSampExpT, firstSampTint=firstSampTint,
                                                 verbose=verbose)
            
            if(is.null(sampProcess)){
              failedSampAnNum.vec<-c(failedSampAnNum.vec,sep.list[[i]]$Analysis[1])
              failedSampID.vec<-c(failedSampID.vec,sep.list[[i]]$Identifier1[1])
              #*
              if(sum(grepl(" ", failedSampID.vec)>0)){
                failedSampID.vec<-gsub(" ","_",failedSampID.vec)
              }
              #*
              failedSampReas.vec<-c(failedSampReas.vec,"NoSamps")
            }else if(length(sampProcess$PeakNr)==0){
              print(paste("No sample peaks detected after processing for Analysis ",sep.list[[i]]$Analysis[1] ))
              # add to remove vector
              failedSampAnNum.vec<-c(failedSampAnNum.vec,sep.list[[i]]$Analysis[1])
              failedSampID.vec<-c(failedSampID.vec,sep.list[[i]]$Identifier1[1])
              #*
              if(sum(grepl(" ", failedSampID.vec)>0)){
                failedSampID.vec<-gsub(" ","_",failedSampID.vec)
              }
              #*
              failedSampReas.vec<-c(failedSampReas.vec,"NoSamps_afterProc")
            }else{
              # check sample iso ratios for d18O/16O
              sampIsoCheck<-isoR_similarityDXF(vend.df=sep.list[[i]],
                                               peakNr.vec=sampProcess$PeakNr,
                                               sdC.thresh=sdCsampIso.thresh,
                                               sdO.thresh=sdOsampIso.thresh,
                                               verbose=verbose)
              
              if((!sampIsoCheck[[1]]$d13C12C)|(!sampIsoCheck[[1]]$d18O16O)){
                sampReasInd.vec<-which(sampIsoCheck[[1]][1,]==FALSE)
                sampReason<-colnames(sampIsoCheck[[1]])[sampReasInd.vec]
                if(length(sampReason)==2){
                  # both isotope SDs failed
                  sampReason<-c("d13C12C_d18O16O")
                }
                failedSampIso.vec<-c(failedSampIso.vec,sep.list[[i]]$Analysis[1])
                failedSampIsoID.vec<-c(failedSampIsoID.vec,sep.list[[i]]$Identifier1[1])
                #*
                if(sum(grepl(" ", failedSampIsoID.vec)>0)){
                  failedSampIsoID.vec<-gsub(" ","_",failedSampIsoID.vec)
                }
                #*
                failedSampIsoSD.vec<-c(failedSampIsoSD.vec,sampIsoCheck[[2]]$`SD_d18O16O`)
                failedSampIsoSDC.vec<-c(failedSampIsoSDC.vec,sampIsoCheck[[2]]$`SD_d13C12C`)
                failedSampIsoReason.vec<-c(failedSampIsoReason.vec,sampReason)
              }else{ # add to sample data if it passes
                samps.list[[sampIndex]]<-sampProcess
                sampIndex<-sampIndex+1
              }
            }
          }
        }
      }
    }
    # return data
    ret.list<-list()
    # create refT check output
    if(is.null(failedRefAnNum.vec)){
      #print("No analyses failed reference times check")
      failedRefAnNum.vec<-c(NA)
      failedRefID.vec<-c(NA)
      failedRefDat.list<-NA
    }
    failedRefDatRetList<-list()
    failedRefAnNum.df<-as.data.frame(matrix(c(failedRefAnNum.vec,failedRefID.vec),ncol=2))
    colnames(failedRefAnNum.df)<-c("failed_refT_Analysis","failed_refT_Identifier1")
    failedRefDatRetList[[1]]<-failedRefAnNum.df
    failedRefDatRetList[[2]]<-failedRefDat.list
    # create refIso check output
    if(is.null(failedRefIso.vec)){
      #print("No analyses failed reference isotope ratio check")
      failedRefIso.vec<-c(NA)
      failedRefIsoID.vec<-c(NA)
      failedRefIsoSD.vec<-c(NA)
      failedRefIsoSDC.vec<-c(NA)
      failedRefIsoReason.vec<-c(NA)
    }
    failedRefIso.df<-as.data.frame(matrix(c(failedRefIso.vec,failedRefIsoID.vec,failedRefIsoSD.vec,failedRefIsoSDC.vec,failedRefIsoReason.vec),ncol=5))
    colnames(failedRefIso.df)<-c("failed_refIso_Analysis","failed_refIso_Identifier1","failed_refIsoSD_d18O16O","failed_refIsoSD_d13C12C","failed_refIsoSD_reason")
    # create refIntensity output
    if(is.null(failedRefInt.vec)){
      #print("No analyses failed reference intensity check")
      failedRefInt.vec<-c(NA)
      failedRefIntID.vec<-c(NA)
      failedRefIntRD.vec<-c(NA)
    }
    failedRefInt.df<-as.data.frame(matrix(c(failedRefInt.vec,failedRefIntID.vec,failedRefIntRD.vec),ncol=3))
    colnames(failedRefInt.df)<-c("failed_refInt_Analysis","failed_refInt_Identifier1","failed_refInt_relDiff")
    # create sample peak processing output
    if(is.null(failedSampAnNum.vec)){
      #print("No analyses failed sample processing")
      failedSampAnNum.vec<-c(NA)
      failedSampID.vec<-c(NA)
      failedSampReas.vec<-c(NA)
    }
    failedSampAnNum.df<-as.data.frame(matrix(c(failedSampAnNum.vec,failedSampID.vec,failedSampReas.vec),ncol=3))
    colnames(failedSampAnNum.df)<-c("failed_sampProc_Analysis","failed_sampProc_Identifier1","failed_SampProc_Reason")
    # create sampIso output
    if(is.null(failedSampIso.vec)){
      #print("No analyses failed sample isotope ratio check")
      failedSampIso.vec<-c(NA)
      failedSampIsoID.vec<-c(NA)
      failedSampIsoSD.vec<-c(NA)
      failedSampIsoSDC.vec<-c(NA)
      failedSampIsoReason.vec<-c(NA)
    }
    failedSampIso.df<-as.data.frame(matrix(c(failedSampIso.vec,failedSampIsoID.vec,failedSampIsoSD.vec,failedSampIsoSDC.vec,failedSampIsoReason.vec),ncol=5))
    colnames(failedSampIso.df)<-c("failed_sampIso_Analysis","failed_sampIso_Identifier1","failed_sampIsoSD_d18O16O","failed_sampIsoSD_d13C12C","failed_sampIsoSD_reason")
    
    # build return data
    ret.list[[1]]<-failedPkNum.df
    ret.list[[2]]<-failedRefDatRetList#failedRefAnNum.df
    ret.list[[3]]<-failedRefIso.df
    ret.list[[4]]<-failedRefInt.df
    ret.list[[5]]<-failedSampAnNum.df
    ret.list[[6]]<-failedSampIso.df
    ret.list[[7]]<-refs.list
    ret.list[[8]]<-samps.list
    return(ret.list)
  }else{
    print("no data passed checks")
    return()
  }
  
}


# (3)
#' writeFailStats: function that writes output from QA/QC to txt files
#' @param processed.list output list from removeFailedAnalysesDXF
#' @param outputID string to identify the current data
#' @param writeDir ; default value = getwd()
#' @examples 
#' Usage example 
#' writeFailStats(processed.list=processed.list,outputID="abiotic",writeDir=getwd())
#' dont export yet
writeFailStats<-function(processed.list,outputID,writeDir=getwd()){
  # analyses that failed the peak number check
  pkNrFileName<-paste(outputID,"_failedPkNr.csv",sep="")
  write.table(processed.list[[1]],pkNrFileName,row.names=F,quote=F,sep=",")
  # analyses that failed the reference peak times check
  refTfileName<-paste(outputID,"_failedRefT.csv",sep="")
  write.table(processed.list[[2]][[1]],refTfileName,row.names=F,quote=F,sep=",")
  # analyses that failed the reference peak isotope ratio test
  refIsoFileName<-paste(outputID,"_failedRefIso.csv",sep="")
  write.table(processed.list[[3]],refIsoFileName,row.names=F,quote=F,sep=",")
  # analyses that failed the reference peak intensity similarity test
  refIntFileName<-paste(outputID,"_failedRefInt.csv",sep="")
  write.table(processed.list[[4]],refIntFileName,row.names=F,quote=F,sep=",")
  # analyses that failed sample peak processing (removal of flush peak and first sample)
  sampProcFileName<-paste(outputID,"_failedSampProc.csv",sep="")
  write.table(processed.list[[5]],sampProcFileName,row.names=F,quote=F,sep=",")
  # analyses that failed sample peak isotope ratio sd test
  sampIsoFileName<-paste(outputID,"_failedSampIso.csv",sep="")
  write.table(processed.list[[6]],sampIsoFileName,row.names=F,quote=F,sep=",")
  if(length(processed.list)==9){
    refSampIntFileName<-paste(outputID,"_failedRefSampInt.csv",sep="")
    write.table(processed.list[[9]],refSampIntFileName,row.names=F,quote=F,sep=",")
  }
}


# (4)
#' separate_by_analysis_numDXF
#' @param vend.df dataframe of the vendor tables from a directory of dxf files
#' @return list of vendor dataframes separated by experiment (analysis number)
#' @examples 
#' Usage example
#' separate_by_analysis_numDXF(vend.df=vend.df)
#' @export
separate_by_analysis_numDXF<-function(vend.df){
  # make lists of data separated by analysis numbers
  analysis_nums<-vend.df$Analysis
  totNumRows<-length(analysis_nums)
  differentAnalyses<-unique(analysis_nums)
  analysisData.list<-list()
  analysisSetInd.vec<-c()
  uniqueAnalysisInd<-1
  # check for end of analysis numbers
  for(i in seq(2,(totNumRows+1))){
    prevAnalysis<-analysis_nums[i-1]
    currAnalysis<-analysis_nums[i]
    # check for end of analysis nums
    if(i==(totNumRows+1)){
      #add last element to vec and then to list
      #print("End of dataset")
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
      analysisData.list[[uniqueAnalysisInd]]<-vend.df[analysisSetInd.vec,]
      break
    }
    # check if the analysis numbers are the same
    if(currAnalysis==prevAnalysis){
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
    } else{ # different analysis numbers
      #print("Reached end of analysis set")
      # add i-1 to vec
      analysisSetInd.vec<-c(analysisSetInd.vec,i-1)
      # add analysis set to list
      analysisData.list[[uniqueAnalysisInd]]<-vend.df[analysisSetInd.vec,]
      # reset analysis index vector and set next list index
      analysisSetInd.vec<-c()
      uniqueAnalysisInd<-uniqueAnalysisInd+1
    }
  }
  return(analysisData.list)
}


# (5)
#' ref_samp_intensity_check
#' don't export yet
ref_samp_intensity_check<-function(currFiltered,
                                   lengthEqual,
                                   dataName="data",
                                   amplName="Ampl44",
                                   relDiffInt.thresh=0.75,
                                   verbose=T){
  if(lengthEqual){
    # separate dfs into lists
    refFiltered.df<-currFiltered[[1]]
    sampFiltered.df<-currFiltered[[2]]
    # separate by experiment
    refFiltered.list<-separate_by_analysis_numDXF(refFiltered.df)
    sampFiltered.list<-separate_by_analysis_numDXF(sampFiltered.df)
    # check peak similarities
    rmAnalysis.vec<-c()
    rmInd.vec<-c()
    
    failedRefSampIntID1.vec<-c()
    failedRefSampRelDiff.vec<-c()
    failedRefSampIntAnalysis.vec<-c()
    for(k in seq(1,length(refFiltered.list))){
      # make vecs of intensities
      analysis<-refFiltered.list[[k]]$Analysis[1]
      ID1<-refFiltered.list[[k]]$Identifier1[1]
      refInt.vec<-refFiltered.list[[k]]$Ampl44 #TODO: make more general
      refPkNr.vec<-refFiltered.list[[k]]$PeakNr
      sampInt.vec<-sampFiltered.list[[k]]$Ampl44 #TODO: make more general
      sampPkNr.vec<-sampFiltered.list[[k]]$PeakNr
      # put data together
      refSampPkNr.vec<-sort(c(refPkNr.vec,sampPkNr.vec))
      refSampVendAmpl.vec<-c()
      refInd=0
      sampInd=0
      for(j in seq(1,length(refSampPkNr.vec))){
        # match peak numbers with intensity amplitudes
        if(refSampPkNr.vec[j] %in% refPkNr.vec){
          refInd<-refInd+1
          refSampVendAmpl.vec<-c(refSampVendAmpl.vec,refInt.vec[refInd])
        }else if(refSampPkNr.vec[j] %in% sampPkNr.vec){
          sampInd<-sampInd+1
          refSampVendAmpl.vec<-c(refSampVendAmpl.vec,sampInt.vec[sampInd])
        }
      }
      # combine in df
      refSampAmpl.mat<-matrix(c(refSampPkNr.vec,refSampVendAmpl.vec),ncol=2)
      refSampAmpl.df<-as.data.frame(refSampAmpl.mat)
      colnames(refSampAmpl.df)<-c("PeakNr",amplName)
      # peak intensity similarity check using intensity_similarityDXF
      refSampIntCheck<-intensity_similarityDXF(vendAmpl=refSampAmpl.df$Ampl44,#TODO: make more general
                                               amplName=amplName,
                                               peakNr.vec=refSampAmpl.df$PeakNr,
                                               relDiffInt.thresh=relDiffInt.thresh)
      
      maxRefSampRelDiff<-refSampIntCheck[[3]]
      # failedRefSampIntID1.vec<-c()
      # failedRefSampRelDiff.vec<-c()
      # failedRefSampIntAnalysis.vec<-c()
      if(!refSampIntCheck[[1]]){
        rmInd=k
        if(verbose==T){print(paste("Significant difference in intensities of reference and sample peaks detected: Analysis ",analysis,sep=""))}
        #print(paste("",maxRefSampRelDiff),sep="")
        # put that analysis number in vec to be removed
        rmAnalysis.vec<-c(rmAnalysis.vec,analysis)
        rmInd.vec<-c(rmInd.vec,rmInd)
        failedRefSampIntAnalysis.vec<-c(failedRefSampIntAnalysis.vec,analysis)
        #failedRefSampIntAnalysis.vec<-rmAnalysis.vec
        failedRefSampIntID1.vec<-c(failedRefSampIntID1.vec,ID1)
        failedRefSampRelDiff.vec<-c(failedRefSampRelDiff.vec,maxRefSampRelDiff)
      }
    }
  }
  # update currFiltered to write new results to file
  if(length(rmAnalysis.vec)>0){
    #remove those analyses from the lists
    #print("removing analyses...")
    refFiltered.list[rmInd.vec]<-NULL
    sampFiltered.list[rmInd.vec]<-NULL
    # change lists back to dfs
    if(length(refFiltered.list)!=0){
      ret.list<-threadRefSampListsToDF(refFiltered.list,sampFiltered.list)
      currFiltered[[1]]<-ret.list[[1]]
      currFiltered[[2]]<-ret.list[[2]]
    }else{ # all data processed out!
      currFiltered[[1]]<-as.data.frame(matrix(NA))
      currFiltered[[2]]<-as.data.frame(matrix(NA))
    }
  }
  
  if(length(failedRefSampIntAnalysis.vec)>0){
    failedRefSampInt.mat<-matrix(c(failedRefSampIntAnalysis.vec,
                                   failedRefSampIntID1.vec,
                                   failedRefSampRelDiff.vec),
                                 ncol=3)
    failedRefSampInt.df<-as.data.frame(failedRefSampInt.mat)
    colnames(failedRefSampInt.df)<-c("failed_RefSampInt_Analysis", "failed_RefSampInt_Identifier1","failed_RefSampInt_RelDiff")
    fileName<-paste(dataName,"_failedRefSampInt.csv",sep="")
    # write failed refSampInt check results to file
    write.table(failedRefSampInt.df,file=fileName,quote=F,row.names=F,sep=",")
  }
  
  return(currFiltered)
}


# (6)
#' intensity_similarityDXF
#' don't export yet
intensity_similarityDXF<-function(vendAmpl,
                                  amplName="Ampl44",
                                  peakNr.vec,
                                  relDiffInt.thresh=0.1,
                                  verbose=T){
  intCheck.list<-list()
  # analyze intensity similarity by relative difference
  int.check<-FALSE
  relDiff<-c()
  rel.vec<-c()
  for(i in seq(2,length(vendAmpl))){
    curr.relDiff<-abs(vendAmpl[i]-vendAmpl[i-1])/max(vendAmpl[i],vendAmpl[i-1])
    rel.vec<-c(rel.vec,curr.relDiff)
    #print(curr.relDiff)
    if(curr.relDiff<relDiffInt.thresh){
      relDiff<-c(relDiff,TRUE)
    }else{
      relDiff<-c(relDiff,FALSE)
    }
  }
  sumRelDiff<-sum(relDiff)
  if(sumRelDiff==length(relDiff)){
    int.check<-TRUE
    #print("intensities within similarity criteria")
    intCheck.list[[1]]<-int.check
  } else{
    int.check<-FALSE
    #if(verbose==T){print("intensities not within similarity criteria")}
    #return(int.check)
    intCheck.list[[1]]<-int.check
  }
  # dataframe of intensity values at peaks
  intCheck.mat<-matrix(c(peakNr.vec,vendAmpl),ncol=2)
  intCheck.df<-as.data.frame(intCheck.mat)
  colnames(intCheck.df)<-c("PeakNr",amplName)
  # add to list
  intCheck.list[[2]]<-intCheck.df
  intCheck.list[[3]]<-max(rel.vec)
  return(intCheck.list)
}


# (7)
#' threadRefSampListsToDF
#' don't export yet
threadRefSampListsToDF<-function(ref.list,sample.list){
  # thread lists into dataframes
  # need total number of elements in the lists
  # sample list
  numSamp<-length(sample.list)
  if(numSamp>0){
    sampSum<-0
    for(i in seq(1,numSamp)){
      sampIndexLength<-dim(sample.list[[i]])[1]
      sampSum<-sampSum+sampIndexLength
    }
    # reference peaks
    numRef<-length(ref.list)
    refSum<-0
    for(i in seq(1,numRef)){
      refIndexLength<-dim(ref.list[[i]])[1]
      refSum<-refSum+refIndexLength
    }
    # df dims: num cols = num colnames; num row=sums
    # initalize references df
    numRefCols<-length(colnames(ref.list[[1]]))
    ref.df<-as.data.frame(matrix(rep(NA,numRefCols*refSum),ncol=numRefCols))
    colnames(ref.df)<-colnames(ref.list[[1]])
    # initalize samples df
    numSampCols<-length(colnames(sample.list[[1]]))
    samp.df<-as.data.frame(matrix(rep(NA,numSampCols*sampSum),ncol=numSampCols))
    colnames(samp.df)<-colnames(sample.list[[1]])
    # loop through lists and build dfs
    # reference df
    numRefs<-length(ref.list)
    refIndex<-1
    for(i in seq(1,numRefs)){
      # grab ref vend
      refDat<-ref.list[[i]]
      numRefPts<-length(refDat$PeakNr)
      # set indices
      refInd.vec<-seq(refIndex,(refIndex+numRefPts-1))
      # add to df
      ref.df[refInd.vec,]<-refDat
      # update ref index
      refIndex<-refInd.vec[length(refInd.vec)]+1
    }
    # sample df
    numSamps<-length(sample.list)
    sampIndex<-1
    for(i in seq(1,numSamps)){
      # grab ref vend
      sampDat<-sample.list[[i]]
      numSampPts<-length(sampDat$PeakNr)
      # set indices
      sampInd.vec<-seq(sampIndex,(sampIndex+numSampPts-1))
      # add to df
      samp.df[sampInd.vec,]<-sampDat
      # update ref index
      sampIndex<-sampInd.vec[length(sampInd.vec)]+1
    }
    refSamp.list<-list()
    refSamp.list[[1]]<-ref.df
    refSamp.list[[2]]<-samp.df
    return(refSamp.list)
  }
}


# (8)
#' reference_times_checkDXF
#' don't export yet
reference_times_checkDXF<-function(vend.df, 
                                   expectedPeak.num=5, 
                                   diff.t=10, 
                                   expectedStart,
                                   expectedRt,
                                   expectedEnd, 
                                   verbose=T){
  curr.analysis<-vend.df$Analysis[1]
  refs.list<-list()
  peak.nums<-as.numeric(vend.df$PeakNr)
  num.peaks<-length(peak.nums)
  # times from vendor table
  start.times<-as.numeric(vend.df$Start)
  start.times
  Rts<-as.numeric(vend.df$Rt)
  end.times<-as.numeric(vend.df$End)
  # intialize loop data
  start.ind<-c()
  rt.ind<-c()
  end.ind<-c()
  ref.ind<-c()
  expectedIndex<-1
  for(i in seq(1,num.peaks)){ # check times and grab reference peaks if detected
    # initialize ref bools
    start.ref<-FALSE
    rt.ref<-FALSE
    end.ref<-FALSE
    # check start
    if(abs(start.times[i]-expectedStart[expectedIndex])<diff.t){ # ??
      # matching start time
      start.ref<-TRUE
      start.ind<-c(start.ind,i)
    }
    # check Rt
    if(abs(Rts[i]-expectedRt[expectedIndex])<diff.t){
      # matching Rt
      rt.ref<-TRUE
      rt.ind<-c(rt.ind,i)
    }
    # check End
    if(abs(end.times[i]-expectedEnd[expectedIndex])<diff.t){
      # matching end time
      end.ref<-TRUE
      end.ind<-c(end.ind,i)
    }
    if(sum(start.ref,rt.ref,end.ref==3)){#all times match for peak as ref
      #print(paste("Start, Rt, and End times match an expected peak for Peak_Nr ",peak.nums[i],sep=""))
      # add to ref peak ind
      ref.ind<-c(ref.ind,i)
      expectedIndex<-expectedIndex+1
    }
  } # end for
  # check if number of reference peaks found matches the expected number
  numRefs<-FALSE
  foundAll<-(length(ref.ind)==expectedPeak.num)
  if(foundAll){   # return reference peak numbers and times
    numRefs<-TRUE
    #print("Expected number of reference peaks detected at expected times.")
    refs.list[[1]]<-numRefs
  } else{
    if(verbose==T){print(paste("Did not detect the expected number of reference peaks at expected times: Analysis ",curr.analysis,sep=""))}
    refs.list[[1]]<-numRefs
  }
  ret<-cbind(vend.df$PeakNr[ref.ind],start.times[ref.ind],Rts[ref.ind],end.times[ref.ind])
  ret.df<-as.data.frame(ret)
  colnames(ret.df)<-c("PkNr","Start","Rt","End")
  refs.list[[2]]<-ret.df
  # vendor data
  refs.list[[3]]<-vend.df[ref.ind,]
  return(refs.list)
}


# (9)
#' maxNumPeaksDXF
#' don't export yet
maxNumPeaksDXF<-function(vend.df,
                         maxExpectedPks=18,
                         expectedNonSampPks=7,
                         verbose=T,
                         extraPeaks=2){
  
  maxRet.list<-list()
  maxRet.df<-as.data.frame(matrix(rep(NA,2),ncol=2))
  colnames(maxRet.df)<-c("NumPks","Analysis")
  num_peaks<-length(vend.df$PeakNr)
  
  extra.vec<-seq(1,extraPeaks)
  acceptedPeakNums.vec<-c(maxExpectedPks-expectedNonSampPks)
  for(i in extra.vec){
    acceptedPeakNums.vec<-c(acceptedPeakNums.vec,maxExpectedPks-expectedNonSampPks-extra.vec[i])
  }
  weirdNum<-!((num_peaks-expectedNonSampPks) %in% acceptedPeakNums.vec)
  # check if data has more than the max expected number of peaks
  #weirdNum<-( (num_peaks-expectedNonSampPks != acceptedPeakNums.vec) && (num_peaks-expectedNonSampPks != 10) && (num_peaks-expectedNonSampPks != 11))
  #if(weirdNum){
  #  print(paste("weird num peaks Analysis ",vend.df$Analysis[1],sep=""))
  #}
  if((num_peaks<=maxExpectedPks) && (num_peaks>0) && (!weirdNum)){
    withinMaxPkNum<-TRUE
    maxRet.list[[1]]<-withinMaxPkNum
    anNum<-vend.df$Analysis[1]
    maxRet.df[1,]<-c(num_peaks,anNum)
    maxRet.list[[2]]<-maxRet.df
    #return(maxRet.list)
  }else{
    if(num_peaks==0){
      print(paste("detected 0 peaks for Analysis ",vend.df$Analysis[1],sep=""))
    }
    withinMaxPkNum<-FALSE
    maxRet.list[[1]]<-withinMaxPkNum
    failedAnNum<-vend.df$Analysis[1]
    if(verbose==T){print(paste("Anomalous number of peaks detected:",num_peaks,"peaks in Analysis",failedAnNum))}
    maxRet.df[1,]<-c(num_peaks,failedAnNum)
    maxRet.list[[2]]<-maxRet.df
  }
  return(maxRet.list)
}


# (10)
#' isoR_similarityDXF
#' don't export yet
isoR_similarityDXF<-function(vend.df,
                             peakNr.vec,
                             sdC.thresh=0.1,
                             sdO.thresh=0.1,
                             verbose=T){
  anaylsis<-vend.df$Analysis[1]
  isoR.list<-list()
  d13C<-vend.df$d13C12C
  d13C.ref<-as.numeric(d13C[peakNr.vec])
  d18O<-vend.df$d18O16O
  d18O.ref<-as.numeric(d18O[peakNr.vec])
  
  # sd of d13C/12C
  sd.13C<-sd(d13C.ref)
  sd.C<-FALSE
  sdnac<-is.na(sd.13C)
  if(!sdnac){
    if(sd.13C<sdC.thresh){
      sd.C<-TRUE
      #print("SD 13C/12C within accepted threshold")
    } else{
      if(verbose==T){print(paste("SD 13C/12C not within accepted threshold: Analysis ",anaylsis,sep=""))}
    }
  } else{
    #print("SD 13C/12C = NA")
  }
  
  # sd of 18O/16O
  sd.18O<-sd(d18O.ref)
  sd.O<-FALSE
  isnasdo<-is.na(sd.18O)
  if(!isnasdo){
    if(sd.18O<sdO.thresh){
      sd.O<-TRUE
      #print("SD 18O/16O within accepted threshold")
    } else{
      if(verbose==T){print(paste("SD 18O/16O not within accepted threshold: Analysis ",anaylsis,sep=""))}
    }
  }else{
    #print("SD 18O/16O = NA")
  }
  
  # dataframe for boolean vals for within threshold
  sdBool.mat<-matrix(c(sd.C,sd.O),ncol=2)
  sdBool.df<-as.data.frame(sdBool.mat)
  colnames(sdBool.df)<-c("d13C12C","d18O16O")
  isoR.list[[1]]<-sdBool.df
  
  # dataframe for numeric value of sds
  sdnum.mat<-matrix(c(sd.13C,sd.18O),ncol=2)
  sdnum.df<-as.data.frame(sdnum.mat)
  colnames(sdnum.df)<-c("SD_d13C12C","SD_d18O16O")
  isoR.list[[2]]<-sdnum.df
  # TODO: use in report
  return(isoR.list)
}


# (11)
#' sample_peaks_processDXF
#' don't export yet
sample_peaks_processDXF<-function(refTimesOutput,
                                  vend.df,
                                  flushExpT=135,
                                  flushTint=15,
                                  firstSampExpT=275,
                                  firstSampTint=15,
                                  verbose=T){
  allPeaks<-seq(1,length(vend.df$PeakNr))
  # start with all as samples then remove reference, flush and 1st sample peak
  samplePeaks<-allPeaks
  # after 4th reference peak
  pkNrFirstSampleAfterRefs<-as.numeric(refTimesOutput[[2]]$PkNr[4])+1
  # grab sample peak vendor info using reference output
  refPeaks<-as.numeric(refTimesOutput[[2]]$PkNr)
  numR<-length(refPeaks)
  # check for samples
  anySamples<-length(refPeaks)<length(allPeaks)
  if(anySamples){ # samples may have flush peak but lack other samples
    for(i in seq(1,numR)){
      if(refPeaks[i] %in% samplePeaks){
        matchInd<-match(refPeaks[i],samplePeaks)
        # remove from samplePeaks
        samplePeaks[matchInd]
        samplePeaks<-samplePeaks[-(matchInd)]
      }
    }
    # get sample peak vendor info
    sample.vend<-vend.df[samplePeaks,]
  }else{
    if(length(refTimesOutput[[3]]$Analysis)>0){
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      if(verbose==T){print(paste("No sample peaks detected for Analysis ",analysisNum,sep=""))}
    }else{
      #if(verbose==T){print("No sample peaks detected")}
    }
    return()
  }
  if(refTimesOutput[[1]]==FALSE){
    # check if analysis number present in vend.df, return if so
    if(length(refTimesOutput[[3]]$Analysis)>0){
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      if(verbose==T){print(paste("FALSE at reference_times_check[[1]]: check reference peak times for quality for Analysis ",analysisNum,sep=""))}
    }else if(length(refTimesOutput[[3]]$Date)<0){
      date<-refTimesOutput[[3]]$Date
      if(verbose==T){print(paste("FALSE at reference_times_check[[1]]: check reference peak times for quality for date ",date,sep=""))}
    } else{
      if(verbose==T){print("FALSE at reference_times_check[[1]]: check reference peak times for quality")}
    }
  }
  # remove flush peak and 1st sample peak
  procSample<-sample.vend
  procSampleStart<-as.numeric(procSample$Start)
  flushInd=0
  
  # find flush peak by time**
  for(i in seq(1,length(procSampleStart))){
    if(abs(flushExpT-procSampleStart[i])<flushTint){ #within 10 secs? or use a different way to ID flush peak?
      flushPeak<-procSample[i,]
      flushInd<-i
    }
  }
  flushNotPresent<-flushInd==0
  if(flushNotPresent){
    #if(verbose==T){print("No flush peak detected within expected time interval")}
  } else{ # remove flush peak
    #print("Flush peak detected within expected time interval")
    procSample<-procSample[-flushInd,]
  }
  # first sample
  firstSampleInd<-which(procSample$PeakNr==pkNrFirstSampleAfterRefs)
  firstSampleVend<-procSample[firstSampleInd,]
  firstSampleTime<-as.numeric(firstSampleVend$Rt)
  # check if other samples present (other than flush)
  if(length(firstSampleInd)==0){
    if(length(refTimesOutput[[3]]$Analysis)>0){ # if there's an analysis number
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      if(verbose==T){print(paste("No sample peaks after flush detected in Analysis ",analysisNum,sep=""))}
    }else{
      #if(verbose==T){print("No sample peaks after flush detected")}
    }
    return()
  }
  # check that the time for the first sample peak is in the expected interval
  sample1withinTime<-abs(firstSampExpT-firstSampleTime)<firstSampTint
  if(sample1withinTime){
    #print("First sample peak detected within expected time interval")
    procSample<-procSample[-firstSampleInd,]
  } else{ # dont remove first sample peak
    analysisNum<-refTimesOutput[[3]]$Analysis[1]
    if(verbose==T){print(paste("First sample peak not within expected time interval: Analysis ",analysisNum,sep=""))}
  }
  # check if there are any sample peaks left after processing out the first 2
  pkNrProc.len<-length(procSample$PeakNr)
  if(pkNrProc.len==0){
    if(length(refTimesOutput[[3]]$Analysis)>0){ # if there's an analysis number
      analysisNum<-refTimesOutput[[3]]$Analysis[1]
      if(verbose==T){print(paste("No sample peaks left after processing in Analysis ",analysisNum,sep=""))}
    }else{
      if(verbose==T){print("No sample peaks left after processing ")}
    }
    return()
  }
  #print(paste("Removed peaks at Rts ",round(as.numeric(flushPeak$Rt),3), "s (flush peak) and ",round(firstSampleTime,3),"s (1st sample peak)",sep=""))
  return(procSample)
}


# (12)
#' all_dxf_files: function that returns the names of all dxf files in the current directory
#' @return vector containing the file names of all dxf files in the current directory
#' @examples 
#' Usage example
#' dxfFileNames <- all_dxf_files()
#' @export
all_dxf_files<-function(){
  all_files<-list.files()
  dxfFiles<-c()
  for(i in seq(1:length(all_files))){
    currFile<-all_files[i] # get file name
    if(grepl(".dxf", currFile)){
      dxfFiles<-c(dxfFiles,currFile)
    }
  }
  return(dxfFiles)
}


# (13)
#' vendor_info: function that returns the vendor table for a given dxf file
#' @param file string that gives the name of the dxf file
#' @return dataframe containing the vendor table informatino 
#' @examples
#' Usage example
#' vend.df <- vendor_info("170526_Na2SO4 L_(1).dxf")
#' @export
vendor_info<-function(file){
  msdat<-iso_read_continuous_flow(file)
  file.info<-msdat %>% iso_get_file_info()
  ident1<-file.info$`Identifier 1`
  vendor_info<-msdat %>% iso_get_vendor_data_table()
  vendor_info.df<-as.data.frame(vendor_info)
  return(vendor_info.df)
}


# (14)
#' vendor_info_all: function that returns the vendor tables for a vector of dxf file names
#' @param files: vector containing dxf file names as elements
#' @return list of vendor table dataframes
#' @examples 
#' Usage example
#' dxf_files.vec <- all_dxf_files()
#' vend.list <- vendor_info_all(dxf_files.vec)
#' @export 
vendor_info_all<-function(files){
  vi.list<-list()
  for(i in seq(1,length(files))){
    vi.df<-vendor_info(files[i])
    vi.list[[i]]<-vi.df
  }
  return(vi.list)
}


# (15)
#' file_info: function that returns the isoreader file_info 
#' @param files 
#' @return file info data frame
#' @examples 
#' Usage example
#' @export
file_info<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  file_info<-msdat %>%
    iso_get_file_info(
      select = c(
        #rename?
        Identifier1 = `Identifier 1`,
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
  file_info.df
  return(file_info.df)
}


# (16)
#' rmRefDatDXF
#' don't export yet
rmRefDatDXF<-function(procList){
  refsFiltered<-procList[[7]]
  sampsFiltered<-procList[[8]]
  analysesEq<-analysisNumsEqualDXF(refsFiltered,sampsFiltered)
  #if(!analysesEq[[1]]){#remove ref data
  rm<-removeRefAnalysisDXF(procList)
  #}
  # check if analysis numbers are all the same
  finalProcRefList<-rm[[1]]
  finalProcSampList<-rm[[2]]
  filter<-threadRefSampListsToDF(finalProcRefList,finalProcSampList)
  return(filter)
}


# (17)
#' analysisNumsEqualDXF
#' don't export yet
analysisNumsEqualDXF<-function(refsList,sampsList){
  lenRefs<-length(refsList)
  if(lenRefs==0){
    print("no data in refs...terminating check")
  }
  else{
    refsAnNums<-c()
    for(i in seq(1,lenRefs)){
      refsAnNums<-c(refsAnNums,as.numeric(refsList[[i]]$Analysis[1]))
    }
    
    lenSamps<-length(sampsList)
    sampsAnNums<-c()
    if(lenSamps>0){
      for(i in seq(1,lenSamps)){
        sampsAnNums<-c(sampsAnNums,as.numeric(sampsList[[i]]$Analysis[1]))
      }
      
      # if they are equal
      lenEq<-(length(refsAnNums)==length(sampsAnNums))
      if(lenEq){
        analysisNumsEq<-sum(refsAnNums-sampsAnNums)==0
        if(analysisNumsEq){
          return(analysisNumsEq)
        }
      }else{ # if they are not equal, return analysis nums that are different?
        anNum.list<-list()
        anRef.df<-as.data.frame(matrix(refsAnNums,ncol=1))
        colnames(anRef.df)<-c("ref_an")
        anS.df<-as.data.frame(matrix(sampsAnNums,ncol=1))
        colnames(anS.df)<-c("samps_an")
        anNum.list[[1]]<-FALSE
        anNum.list[[2]]<-anRef.df
        anNum.list[[3]]<-anS.df
        return(anNum.list)
      }
    }
  }
}


# (18)
#' removeRefAnalysisDXF
#' don't expot yet
removeRefAnalysisDXF<-function(filtered.list, refInd=7, sampInd=8,
                               sampProcFailedInd=5, sampIsoFailInd=6){
  retSR.list<-list()
  refs<-filtered.list[[refInd]]
  samps<-filtered.list[[sampInd]]
  
  failedSampIsoAnNum.vec<-filtered.list[[sampIsoFailInd]][,1]
  failedSampProcAnNum.vec<-filtered.list[[sampProcFailedInd]][,1]
  failedSampAn.vec<-c(failedSampIsoAnNum.vec,failedSampProcAnNum.vec)
  
  removeInd<-c()
  lenRefsList<-length(refs)
  if(lenRefsList>0){
    for(i in seq(1,lenRefsList)){
      anNum<-refs[[i]]$Analysis[1]
      if(anNum %in% failedSampAn.vec){
        removeInd<-c(removeInd,i)
      }
    }
    ## remove all at once
    refs[removeInd]<-NULL #closes the hole
    
    retSR.list[[1]]<-refs
    retSR.list[[2]]<-samps
    
    return(retSR.list)
  }
}



# (19)
#' QAQC_IRMS: wrapper function for QA/QC that writes and returns quality IRMS data and QC check diagnostics
#' @param unfilteredPath path to dxf files for QA/QC; default = current working directory
#' @param expRef.df dataframe of expected reference peak times with column names = c("Ref_Peak_Nr", "Expected_Start", "Expected_Rt", "Expected_End")
#' @param diff.t time interval in seconds for reference peak retention times; default = 10
#' @param checkIntStand whether to perform internal standards calibration; default = F; not currently functional
#' @param internalStandID vector of internal standard names for calibration; default = c("L1","H1","LW")
#' @param useColNames vector of column names to extract vendor data from in the dxf files; default = c("fileId","Identifier1","Analysis","Preparation","DateTime",
#'                                   "PeakNr","Start","Rt","End","Ampl44","Ampl45",
#'                                   "Ampl46","BGD44","BGD45","BGD46","rIntensity44","rIntensity45",
#'                                   "rIntensity46","rIntensityAll","Intensity44","Intensity45",
#'                                   "Intensity46","IntensityAll","ListFirstPeak","rR45CO244CO2",
#'                                   "rR46CO244CO2","IsRef","R45CO244CO2","RefName","rd45CO244CO2",
#'                                   "d45CO244CO2", "R46CO244CO2", "rd46CO244CO2","d46CO244CO2",
#'                                   "R13C12C","d13C12C","AT13C12C","R18O16O","d18O16O", "AT18O16O",
#'                                   "R17O16O","d17O16O","Rps45CO244CO2","Rps46CO244CO2"),
#'                                  outputID = "abiotic")
#' @param dataName string giving the name of the current dataset
#' @param maxPkNum maximum number of accepted peaks in an experiment; default = 18 
#' @param expectedNonSampPks expected number of non-sample peaks in an experiment; default = 7
#' @param sdCrefIso.thresh maximum acceptable delta 13C for the reference peaks; default value = 0.1
#' @param sdOrefIso.thresh maximum acceptable delta 18O for the reference peaks; default value = 0.1
#' @param checkRelDiffIntensity whether to check the reference:sample peak relative difference in peak intensities; default = T
#' @param amplName name of the column to use for the peak intensity checks; default value = "Ampl44"
#' @param relDiffInt.thresh maximum acceptable difference in reference peak relative intensity; default value = 0.1
#' @param sdCsampIso.thresh maximum acceptable delta 13C for the sample peaks; default value = 0.3
#' @param sdOsampIso.thresh maximum acceptable delta 18O for the sample peaks; default value = 0.2
#' @param verbose whether to print information about the checks; default value = T
#' @return list of two elements: ret.list[[1]] - reference peak data that passed QA/QC
#'                               ret.list[[2]] - sample peak data that passed QA/QC
#' @examples 
#' Usage example
#' qc_data.list <- QAQC_IRMS()
#' @export 
QAQC_IRMS<-function(unfilteredPath, 
                    expRef.df, 
                    diff.t=10,
                    checkIntStand=F, 
                    standAcceptedVals.vec=c(-8.55,4.85,-3.85),
                    standAcceptedSD.vec=c(0.2,0.2,0.2),
                    internalStandID=c("L1","H1","LW"),
                    useColNames=c("fileId","Identifier1","Analysis","Preparation","DateTime",
                                  "PeakNr","Start","Rt","End","Ampl44","Ampl45",
                                  "Ampl46","BGD44","BGD45","BGD46","rIntensity44","rIntensity45",
                                  "rIntensity46","rIntensityAll","Intensity44","Intensity45",
                                  "Intensity46","IntensityAll","ListFirstPeak","rR45CO244CO2",
                                  "rR46CO244CO2","IsRef","R45CO244CO2","RefName","rd45CO244CO2",
                                  "d45CO244CO2", "R46CO244CO2", "rd46CO244CO2","d46CO244CO2",
                                  "R13C12C","d13C12C","AT13C12C","R18O16O","d18O16O", "AT18O16O",
                                  "R17O16O","d17O16O","Rps45CO244CO2","Rps46CO244CO2"), 
                    dataName,
                    maxPkNum=18, 
                    expectedNonSampPks=7,
                    sdCrefIso.thresh=0.1,
                    sdOrefIso.thresh=0.1,
                    checkRelDiffIntensity=T,
                    amplName="Ampl44",
                    relDiffInt.thresh=0.1,
                    sdCsampIso.thresh=0.3,
                    sdOsampIso.thresh=0.2,
                    verbose=T #TODO: save "bad" data
){
  # return data that passes QC
  ret.list<-list()
  
  ### QC for DXF files
  # create directory for results
  writeDir<-paste(unfilteredPath,"QCresults/",sep="")
  
  combList<-combineVendFileInfo(unfilteredPath,useColNames,dataName)
  origNumFiles<-length(combList) 
  
  if(verbose==T){
    print(paste("original number of files: ",length(combList)))
  }
  
  # filter out bad dxf data
  currProc<-removeFailedAnalysesDXF(sepList=combList,
                                    expRef.df=expRef.df,
                                    maxPkNum=maxPkNum,
                                    expRefPkNr=dim(expRef.df)[1],
                                    diff.t=diff.t,
                                    sdCrefIso.thresh=sdCrefIso.thresh, 
                                    expectedNonSampPks=expectedNonSampPks,
                                    relDiffInt.thresh=0.1, amplName=amplName,
                                    sdOrefIso.thresh=sdOrefIso.thresh,
                                    sdCsampIso.thresh=sdCsampIso.thresh,
                                    sdOsampIso.thresh=sdOsampIso.thresh,
                                    verbose=verbose)
  
  # adjust reference data for any samples removed
  # list of two lists: one for ref peak dfs, one for samps
  currFiltered<-rmRefDatDXF(currProc)
  # want same number of analyses in both
  lengthEqual<-length(currFiltered[[1]])==length(currFiltered[[2]])
  #lengthEqual
  if(length(currFiltered)>0){
    # check relative differences bt samp and ref intensity
    if(checkRelDiffIntensity && lengthEqual){
      #print("checking relative difference between sample and reference peak intensity")
      # TODO: 
      #currFiltered,lengthEqual,dataName="data",amplName="Ampl44",relDiffInt.thresh=0.75,verbose=T
      currFiltered<-ref_samp_intensity_check(currFiltered=currFiltered,
                                             lengthEqual=lengthEqual,
                                             dataName=dataName,
                                             amplName=amplName,
                                             relDiffInt.thresh=0.75, #TODO: make another var for this func
                                             verbose=verbose)
    }
  }
  # write to file in curr dir and combined directory
  # reference data
  write.csv(currFiltered[[1]],
            file=paste(dataName,"_refs_allChecks.csv",sep=""),quote=F,row.names=F)
  # sample data
  write.csv(currFiltered[[2]],file=paste(dataName,"_samps_allChecks.csv",sep=""),quote=F,row.names=F)
  # write failure summary
  stat<-writeFailStats(currProc,dataName,writeDir)
  
  sampsFiltered.df<-currFiltered[[2]]
  
  ### internal standards check: optional flag
  if(checkIntStand==T){
    print("analyzing internal standards...")
    # get avg/sd of internal standards, do linear fit for calibration 
    int_stand.list <- internal_standards_summary(data.df=sampsFiltered.df, dataName=dataName,
                                                 standName.vec=internalStandID,
                                                 standAcceptedVals.vec=standAcceptedVals.vec,
                                                 standAcceptedSD.vec=standAcceptedSD.vec)
    currFiltered[[3]]<-int_stand.list
    
  }
 
  print(paste("original number of experiments: ",
              length(combList)),sep="")
  print(paste("QA/QC'd number of experiments: ",
              length(unique(currFiltered[[2]]$Analysis))),sep="")
  return(currFiltered)
  
}


# (20)
#' read_summary: read and print a summary of mass spec data from a .dxf file
#' @param filename character string of the name of the .dxf file of data
#' @return summary table of file contents
#' @examples
#' Usage example
#' read_summ("170525_NaHCO3 L + NaCl L_.dxf")
#' @export
read_summary<-function(filename){
  msdat<-iso_read_continuous_flow(filename)
  summ<-summary(msdat)
  print(summ)
}


# (21)
#'raw_data: get the raw data from the specified .dxf file as a dataframe
#' @param file character string defining the file to extract raw data from
#' @return dataframe containing the raw data in the .dxf file
#' @examples
#' Usage example
#' raw_data(data_files)
#' @export
raw_data<-function(file){
  #num_files<-length(files)
  msdat<-iso_read_continuous_flow(file)#files[1:num_files]
  raw_dat<- msdat %>% iso_get_raw_data()
  raw_dat.df<-as.data.frame(raw_dat)
}


# (22)
#' raw_data_all: function to get all raw data from multiple files
#' @param files vector of filenames to get raw data from
#' @return a list containing raw data for each file
#' @examples
#' Usage example
#' rawList<-raw_data_all(files)
#' @export
raw_data_all<-function(files){
  raw.list<-list()
  for(i in seq(1,length(files))){
    raw.dat<-raw_data(files[i])
    raw.list[[i]]<-raw.dat
  }
  return(raw.list)
}


# (23)
#' resistor_data: get resistor info for .dxf files as a dataframe
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of resistor info in the .dxf files
#' @examples
#' Usage example
#' resistor_df(data_files)
#' @export
resistor_data<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  resistors<-iso_get_resistors(msdat)
  resistors.df<-as.data.frame(resistors)
}


# (24)
#' reference_values_ratio: get isotopic reference values including ratios
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of reference values in the .dxf files
#' @examples
#' Usage example
#' reference_values_ratio(data_files)
#' @export
reference_values_ratio <- function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  # reference values with ratios
  ref_with_ratios<-msdat %>% iso_get_standards()
  ref_with_ratios.df<-as.data.frame(ref_with_ratios)
}


# (25)
#' reference_values_no_ratio: get isotopic reference values without ratios
#' @param files vector containing character strings of .dxf file names
#' @return dataframe of reference values in the .dxf files
#' @examples
#' Usage example
#' reference_values_no_ratio(data_files)
#' @export
reference_values_no_ratio<-function(files){
  num_files<-length(files)
  msdat<-iso_read_continuous_flow(files[1:num_files])
  # reference delta values without ratio values
  delta_no_ratio<-msdat %>% iso_get_standards(file_id:reference)
  delta_no_ratio.df<-as.data.frame(delta_no_ratio)
}


# (26)
#' trap_area_allPks: function to calculate area of all peaks in an experiment via trapz (wrapper function for all_PA_trap())
#' @param raw.df dataframe of raw data
#' @param vend.df dataframe of vendor data
#' @param mV.rawName character string of the name of intensity column to use (in raw.df)
#' @return dataframe containing peak numbers and areas
#' @examples
#' Usage Example
#' pkAreas.df<-trap_area_allPks(raw.df=raw.df,vend.df=vi.df,mV.rawName="v44.mV")
#' @export
trap_area_allPks<-function(raw.df,vend.df,mV.rawName){
  # get mass voltage index
  massInd<-which(colnames(raw.df)==mV.rawName)
  massV<-raw.df[,massInd]
  massT<-raw.df$time.s
  massSt<-as.numeric(vend.df$Start)
  massEnd<-as.numeric(vend.df$End)
  massPkNr<-as.numeric(vend.df$Nr.)
  
  rawMass.mat<-matrix(c(massT,massV),ncol=2)
  rawMass.df<-as.data.frame(rawMass.mat)
  colnames(rawMass.df)<-c("time.s",mV.rawName)
  
  vendMass.mat<-matrix(c(massSt,massEnd,massPkNr),ncol=3)
  vendMass.df<-as.data.frame(vendMass.mat)
  colnames(vendMass.df)<-c("Start","End","Pk_Nr")
  
  ret.list<-list()
  ret.list[[1]]<-rawMass.df
  ret.list[[2]]<-vendMass.df
  
  rawIntTime<-ret.list[[1]]
  vendTimePk<-ret.list[[2]]
  
  start<-vendTimePk$Start
  end<-vendTimePk$End
  PkNr<-vendTimePk$Pk_Nr
  time<-rawIntTime$time.s
  mV<-rawIntTime[,2]
  
  # call area func
  allA<-all_PA_trap(start,end,time,mV,PkNr)
  return(allA)
}


# (27)
#' all_PA_trap: function to calculate all peak areas using the integration trapezoidal rule (via trapz)
#' @param start.vec numeric vector of peak start times
#' @param end.vec numeric vector of peak end times
#' @param time.vec numeric vector of times from the raw data
#' @param int.vec numeric vector of intensities from the raw data
#' @return a numeric vector containing all the peak areas in V*s
#' @examples
#' Usage Example
#' all_PA_trap(start.v1,end.v1,time.s,v44,pk.Nrs)
#' dont export yet
all_PA_trap<-function(start.vec,end.vec,time.vec,int.vec,pk.Nrs){
  all_areas<-c()
  for(i in seq(1:length(start.vec))){
    all_areas<-c(all_areas,peak_area_trap(start.vec[i],end.vec[i],time.vec,int.vec))
  }
  all_areas_Vs<-all_areas
  areas.mat<-matrix(c(pk.Nrs,all_areas_Vs),ncol=2)
  areas.df<-as.data.frame(areas.mat)
  colnames(areas.df)<-c("Pk_Nr","trap_area")
  return(areas.df)
}


# (28)
#' peak_area_trap: function that returns the area of a peak using the integration trapezoidal rule via trapz
#' @param start.t start time of the peak
#' @param end.t end time of the peak
#' @param time.vec time vector of times in the interval for the peak
#' @param int.vec intensity vector of intensity values for the peak
#' @return peak area in V*s
#' @examples
#' Usage example
#' peak_area_trap(start.v1[1],end.v1[1],time.s,v45)
#' @export
peak_area_trap<-function(start.t,end.t,time.vec,int.vec){
  # get times for peak
  peak.t<-c()
  time.ind<-c()
  # get peak times and indices
  for(i in seq(1:length(time.vec))){
    if((time.vec[i]>=start.t) && (time.vec[i]<=end.t)){
      peak.t<-c(peak.t,time.vec[i])
      time.ind<-c(time.ind,i)
    }
  }
  # get peak intensities
  peak.mv<-c()
  for(i in seq(1:length(time.ind))){
    peak.mv<-c(peak.mv,int.vec[time.ind[i]])
  }
  #plot(peak.t,peak.mv,type="l")
  peak.df<-as.data.frame(cbind(peak.t,peak.mv))
  colnames(peak.df)<-c("time","intensity")
  #peak.area<-cha(peak.t,peak.mv)
  peak.area<-trapz(peak.t,peak.mv)/1000
  return(peak.area)
}


# (29)
#' generic_plot_all_raw: function that plots all the raw data in a list containing raw data from different experiments
#' @param raw.list list whose elements are the raw.df to be plotted
#' @examples
#' Usage Example
#' generic_plot_all_raw(rawList)
#' @export
generic_plot_all_raw<-function(raw.list){
  raw.length<-length(raw.list)
  par(mfrow=c(2,3))
  for(i in seq(1,raw.length)){
    # get title
    raw.dat<-raw.list[[i]]
    raw.title<-raw.dat$file_id[1]
    raw.title<-paste("Index: ",i," Raw Data ",raw.title,sep="")
    raw.title
    # plot
    if(dim(raw.dat)[1]>0){
      generic_raw_plot(raw.dat,raw.title)
    }
  }
  par(mfrow=c(1,1))
}


# (30)
#' generic_raw_plot: function to plot raw data using generic plot
#' @param raw.df full dataframe of raw data
#' @param title title of the plot
#' @examples
#' Usage Example
#' generic_raw_plot(rawDat.df,"170525_NaHCO3 L + NaCl L_.dxf")
#' @export
generic_raw_plot<-function(raw.df,title){
  #par(mfrow=c(1,1))
  # get intensity data for each mass
  v44<-raw.df$v44.mV
  v45<-raw.df$v45.mV
  v46<-raw.df$v46.mV
  # get time data
  time.s<-raw.df$time.s
  # this way, have to plot the tallest one first to get all in same window
  plot(time.s,v46,type="l",col="magenta",ylab="Intensity (mV)",xlab="Time(s)")#v44_v45_v46.mV
  lines(time.s,v45,type="l",col="green")
  lines(time.s,v44,type="l",col="blue")
  title(main=title)
}


# (31)
#' DXFvendListFromID1: function
#' don't export yet
DXFvendListFromID1<-function(all.df,standName.vec=c("L1","H1","LW")){
  vendRet.list<-list()
  listIndex<-1
  # separate into different experiments using analysis numbers
  sepList<-separate_by_analysis_numDXF(all.df)
  # look just for identifier1 that are in standName.vec
  lenList<-length(sepList)
  for(i in seq(1,lenList)){
    currID<-sepList[[i]]$Identifier1[1]
    currID
    foundStand<-currID %in% standName.vec
    foundStand
    if(foundStand){ # get vend data and add to output list: Analysis,Identifier1,PkNr,d18O16O,d13C/12C
      standVend<-sepList[[i]]
      currAn.vec<-standVend$Analysis
      currID.vec<-standVend$Identifier1
      currPkNr.vec<-standVend$PeakNr
      currd18O<-standVend$d18O16O
      currd13C<-standVend$d13C12C
      curr.df<-as.data.frame(matrix(c(currAn.vec,currID.vec,currPkNr.vec,currd18O,currd13C),ncol=5))
      colnames(curr.df)<-c("Analysis","Identifier1","PkNr","d18O16O","d13C12C")
      # add to list
      vendRet.list[[listIndex]]<-curr.df
      listIndex<-listIndex+1
    }
  }
  return(vendRet.list)
}


# (32)
#' avg_sd_d18O_standards: function to find d18O/16O averages and standard deviations for each set of files represented in input list
#' @param allStandards_d18O.list list built from output of d18O_samples_list(), which retrieves sample peaks and d18O/16O data from specified files
#' @param standNames character vector containing the names of the standards
#' @param standAcceptedVals.vec numeric vector containing the accepted values of the standards
#' @param accStandRatioSD numeric vector containing the accepted standard deviation values for each of the standards (in the same order as names in standNames)
#' @return list of length 3 with the following elements:
#'          [[1]]: dataframes of d18O/16O averages and standard deviations for each set of files from input
#'          [[2]]: accepted and measured d18O/16O dataframe
#'          [[3]]: acceptable and calculated standard deviations for d18O/16O
#' @examples
#' Usage Example
#' L1_d18O.list<-d18O_samples_list(L1fileNames.vec)
#' H1_d18O.list<-d18O_samples_list(H1fileNames.vec)
#' LW_d18O.list<-d18O_samples_list(LWfileNames.vec)
#' allStand_d18O.list<-list()
#' allStand_d18O.list[[1]]<-L1_d18O.list
#' allStand_d18O.list[[2]]<-H1_d18O.list
#' allStand_d18O.list[[3]]<-LW_d18O.list
#' AvgSD_d18O<-avg_sd_d18O_standards(allStand_d18O.list)
#' dont export yet
avg_sd_d18O_standards<-function(allStandards_d18O.list,
                                standNames=c("L1","H1","LW"), isoName="d18O16O",
                                standAcceptedVals.vec=c(-8.55,4.85,-3.85),
                                accStandRatioSD=c(0.2,0.2,0.2)){
  
  # make df for avgs (avgs of all d18O in all files for particular standard)
  standAccepted.mat<-matrix(standAcceptedVals.vec,nrow=length(standAcceptedVals.vec))
  standAccepted.df<-as.data.frame(standAccepted.mat)
  dim(standAccepted.df)
  # start df -- will add average of all measured values
  d18OstandVals_Avgs.df<-data.frame(matrix(rep(NA,3*(length(standNames))),ncol=3))
  colNames<-c(paste("accepted_",isoName,sep=""),paste("measured_",isoName,sep=""))
  colnames(d18OstandVals_Avgs.df)<-c("standard",colNames)
  #rownames(d18OstandVals_Avgs.df)<-standNames
  d18OstandVals_Avgs.df[,1]<-standNames
  d18OstandVals_Avgs.df[,2]<-standAccepted.df
  
  # make df for SDs
  SDAccepted.mat<-matrix(accStandRatioSD,ncol=1)
  SDAccepted.df<-as.data.frame(SDAccepted.mat)
  dim(SDAccepted.df)
  # start df -- will add SD of all measured values
  d18OstandVals_SDs.df<-data.frame(matrix(rep(NA,3*(length(standNames))),ncol=3))
  sdColNames<-c(paste("accepted_SD_",isoName,sep=""),paste("calculated_SD_",isoName,sep=""))
  colnames(d18OstandVals_SDs.df)<-c("standard",sdColNames)
  #rownames(d18OstandVals_SDs.df)<-standNames
  d18OstandVals_SDs.df[,1]<-standNames
  d18OstandVals_SDs.df[,2]<-accStandRatioSD
  
  # initialize lists for return vals
  #ret_files_d18.list<-list()
  ret.list<-list()
  
  # intStandsTotal.list<-list()
  # for(i in seq(1,length(standNames))){
  #   intStandsTotal.list[[i]]<-NA
  # }
  
  d18O.list<-allStandards_d18O.list
  list.len<-length(d18O.list)
  
  # for each file - do at end?
  #alld18OFile.list<-list()
  #alld18OFile.index<-1
  d18OFile.mat<-matrix(rep(NA,list.len*4),ncol=4) 
  d18OFile.df<-as.data.frame(d18OFile.mat)
  fileColNames<-c(paste("avg_",isoName,sep=""),paste("sd_",isoName,sep=""))
  colnames(d18OFile.df)<-c("Analysis","standard",fileColNames)
  
  for(i in seq(1,list.len)){
      # get sample peaks d18O/16O for one file
      dlist.df<-d18O.list[[i]]
      # get the filename
      #fileN<-dlist.df$file_id[1]
      curr_analysis<-dlist.df$Analysis[1]
      curr_stand<-dlist.df$Identifier1[1]
      # get the d18O data for the file
      ind <- which(colnames(dlist.df)==isoName)
      d18O.vec<-as.numeric(dlist.df[,ind])
      #d13C.vec<-as.numeric(dlist.df$d13C12C)
      
      # put in vector of all d18O for the given standard (i iteration, over the 3 files for a standard)
      #alld18O.vec<-c(alld18O.vec,d18O.vec)
      # get avg d18O for each file
      avg<-mean(d18O.vec)
      # get sd for each file
      sd<-sd(d18O.vec)
      # put file_id, avg, and sd for each file in df
      d18OFile.df[i,]<-cbind.data.frame(curr_analysis,curr_stand,avg,sd)
  }
  
  stand_dir.df<-as.data.frame(matrix(rep(NA,3*length(standNames)),ncol=3))
  colnames(stand_dir.df)<-c("standard",fileColNames)
  # for each standard, total avg/sd in directory 
  for(j in seq(1,length(standNames))){
    curr.name<-standNames[j]
    curr.ind<-which(d18OFile.df$standard == curr.name)
    curr.dat <- d18OFile.df[curr.ind,]
    curr.avg <- mean(curr.dat[,3])
    curr.sd <- mean(curr.dat[,4])
    curr.df<-cbind.data.frame(curr.name,curr.avg,curr.sd)
    stand_dir.df[j,]<-curr.df
  }
  stand_dir.df
  
  d18OstandVals_SDs.df[,3] <- stand_dir.df[,3]
  d18OstandVals_Avgs.df[,3] <- stand_dir.df[,2]
  
  ret.list[[1]] <- d18OFile.df
  ret.list[[2]] <- d18OstandVals_Avgs.df
  ret.list[[3]] <- d18OstandVals_SDs.df
  
  # add to list for dfs of file sets
    #alld18OFile.list[[alld18OFile.index]]<-d18OFile.df
    #alld18OFile.index<-alld18OFile.index+1
    
    # avg of all d18O in files for given standard
    #d18OstandAvg<-mean(alld18O.vec)
    # put in df to return
    #d18OstandVals_Avgs.df$`measured_d18O16O`[i]<-d18OstandAvg
    
    # SD
    #sdAllStandFiles<-sd(alld18O.vec)
    #d18OstandVals_SDs.df$`calculated_SD_d18O16O`[i]<-sdAllStandFiles
    
    # return d18O.df in list
    #ret_files_d18.list[[i]]<-d18OFile.df

  # add values to return list
  #ret_d18.list[[1]]<-ret_files_d18.list
  #ret_d18.list[[2]]<-d18OstandVals_Avgs.df # all avgs
  #ret_d18.list[[3]]<-d18OstandVals_SDs.df # all SDs
  
  return(ret.list)
}


# (33)
#' stand_lm: function to perform the linear regression for the internal standards quality check and plot the results
#' @param acceptedMeas.df dataframe with accepted and measured values for internal standards (rownames are thestandard names, colnames=c(accepted,measured))
#' @return list whose first element is the lm model and the second is the model summary
#' @examples
#' Usage Example
#' avgSD_d18O<-avg_sd_d18O_standards(allStand_d18O.list)
#' standDat<-testAvg_d18O[[2]]
#' stand.lm<-stand_lm(standDat)
#' dont export yet
stand_lm<-function(acceptedMeas.df,dataName){
  ret.list<-list()
  # linear regression
  standard.fit<-lm(accepted_d18O16O~measured_d18O16O,data=acceptedMeas.df)
  standard.fit
  fit.summ<-summary(standard.fit)
  # plot
  plot(acceptedMeas.df$measured,acceptedMeas.df$accepted,col=c("orange","green","blue"),
       main=paste("lm for accepted vs measured d18O standards: ",dataName,"\n",
                  "r^2 = ",round(fit.summ$r.squared,4),sep=""), 
       xlab="measured",ylab="accepted",pch=20)
  legend("bottomright", legend=c("L1","LW","H1"),
        col=c("orange","blue","green"), pch=20)
  abline(standard.fit,col="light blue")
  
  ret.list[[1]]<-standard.fit
  ret.list[[2]]<-fit.summ
  return(ret.list)
}


# (34)
#' internal_standards_summary: function 
#' @param data.df dataframe of sample peak data 
#' @param dataName dataset name
#' @param standName.vec vector of internal standard names
#' @param standAcceptedVals.vec vector of internal standard delta values
#' @param standAcceptedSD.vec vector of acceptable internal standard SD in delta values
#' @return list of four elements:
#'   ret.list[[1]]: list of deltas by experiment
#'   ret.list[[2]]: average values of deltas in the dataset
#'   ret.list[[3]]: mean SD of deltas in the dataset
#'   ret.list[[4]]: results of linear model using accepted and measured values
#' Usage example
#' @export
internal_standards_summary <- function(data.df, dataName,
                                       standName.vec=c("L1","H1","LW"),
                                       standAcceptedVals.vec=c(-8.55,4.85,-3.85),
                                       standAcceptedSD.vec=c(0.2,0.2,0.2)){
  ret.list<-list()
  
  standIsoR.list<-DXFvendListFromID1(all.df=samps.dat,standName.vec=standName.vec)
  
  standAvgSD<-avg_sd_d18O_standards(allStandards_d18O.list = standIsoR.list,
                                    standNames=standName.vec,
                                    standAcceptedVals.vec=standAcceptedVals.vec,
                                    accStandRatioSD=standAcceptedSD.vec)
  # first element: summary of avg d18O/16O and SD d18O/16O for all internal standards
  standAnalyses.df<-standAvgSD[[1]]
  write.table(standAnalyses.df,"interal_standards_analysis.csv",row.names=F,quote=F,)   
  # average
  avg_d18O<-standAvgSD[[2]]
  write.table(avg_d18O,"avg_interal_standards.csv",row.names=F,quote=F,sep=",")
  # standard deviations
  sd_d18O<-standAvgSD[[3]]
  write.table(sd_d18O,"sd_internal_standards.csv",row.names=F,quote=F,sep=",")
  
  ## save graph to file
  pdf(file=paste("int_stand_lm.pdf",sep=""),width=6,height=4)
  standlm<-stand_lm(avg_d18O,dataName=dataName)
  dev.off()
  
  ## write stats
  # slope, intercept, r^2
  lmr2<-standlm[[2]]$r.squared
  lmresid<-standlm[[2]]$residuals
  
  lmCoeff.df<-as.data.frame(matrix(c(standlm[[1]]$coefficients,lmr2,lmresid),ncol=6))
  colnames(lmCoeff.df)<-c("intercept","slope","r^2","L1_resid","H1_resid","LW_resid")
  
  # write lm stats to file
  write.table(lmCoeff.df,"intStand_lm.csv",row.names=F,quote=F,sep=",")
  
  # std error, t val and p val
  lmCoeffStat<-standlm[[2]]$coefficients
  write.table(lmCoeffStat,"intStand_coefficients.csv",row.names=T,quote=F,sep=",")
  
  # return data
  ret.list[[1]]<-standIsoR.list
  ret.list[[2]]<-avg_d18O
  ret.list[[3]]<-sd_d18O
  ret.list[[4]]<-standlm
  return(ret.list)
}
