### functions to create labeled ML datasets from QC'd IRMS data

# replace space prep
replace_space_prep<-function(data.df){
  prep.vec<-data.df$Preparation
  newPrep.vec<-c()
  for(i in seq(1,length(prep.vec))){
    currPrep<-prep.vec[i] # get file name
    if(grepl(" ",currPrep)){
      newPrep<-gsub(" ","_",currPrep) # substitute spaces with underscores
      newPrep.vec<-c(newPrep.vec,newPrep)
    }else{
      newPrep.vec<-c(newPrep.vec,currPrep)
    }
    if((grepl("He-Only",currPrep))){ # "He only", "He-Only"
      newPrep<-"He_only"
      newPrep.vec[i]<-newPrep
    }
  }
  data.df$Preparation<-newPrep.vec
  return(data.df)
}


# replace space Identifier1
replace_space_Identifier1<-function(data.df){
  id1.vec<-data.df$Identifier1
  newID1.vec<-c()
  for(i in seq(1,length(id1.vec))){
    currID<-id1.vec[i] # get file name
    if(grepl(" ",currID)){
      newID1<-gsub(" ","_",currID) # substitute spaces with underscores
      newID1.vec<-c(newID1.vec,newID1)
    }
    else{
      newID1.vec<-c(newID1.vec,currID)
    }
  }
  data.df$Identifier1<-newID1.vec
  return(data.df)
}

# calculate peak areas using trapezoid rule
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
  #head(rawMass.df)
  vendMass.mat<-matrix(c(massSt,massEnd,massPkNr),ncol=3)
  vendMass.df<-as.data.frame(vendMass.mat)
  colnames(vendMass.df)<-c("Start","End","PkNr")

  ret.list<-list()
  ret.list[[1]]<-rawMass.df
  ret.list[[2]]<-vendMass.df

  rawIntTime<-ret.list[[1]]
  vendTimePk<-ret.list[[2]]

  start<-vendTimePk$Start
  end<-vendTimePk$End
  PkNr<-vendTimePk$PkNr
  time<-rawIntTime$time.s
  mV<-rawIntTime[,2]

  # call area func
  allA<-all_PA_trap(start,end,time,mV,PkNr)
  return(allA)
}

# dependency func for trap_area_allPks
all_PA_trap<-function(start.vec,end.vec,time.vec,int.vec,pk.Nrs){
  all_areas<-c()
  for(i in seq(1:length(start.vec))){
    all_areas<-c(all_areas,peak_area_trap(start.vec[i],end.vec[i],time.vec,int.vec))
  }
  all_areas_Vs<-all_areas
  areas.mat<-matrix(c(pk.Nrs,all_areas_Vs),ncol=2)
  areas.df<-as.data.frame(areas.mat)
  colnames(areas.df)<-c("PkNr","trapArea")
  return(areas.df)
}


# # add d18O/13C label
d18O13Clab<-function(data.df){
  dOC.vec<-c()
  for(i in seq(1,dim(data.df)[1])){
    d18O13C<-data.df$d18O16O[i]/data.df$d13C12C[i]
    dOC.vec<-c(dOC.vec,d18O13C)
  }
  # add to df
  data.df$d18O13C<-dOC.vec
  return(data.df)
}


# find the averages of the deltas
avg3Deltas<-function(data.df){
  # separate data into list
  data.list<-separate_by_analysis(data.df)
  avgd18O16O.vec<-c()
  avgd13C12C.vec<-c()
  avgd18O13C.vec<-c()
  analysis.vec<-c()
  for(i in seq(1,length(data.list))){ #231 for dxf data
    curr.df<-data.list[[i]]
    analysis.vec<-c(analysis.vec,curr.df$Analysis[1])
    # averages
    avgd18O16O<-mean(curr.df$d18O16O)
    avgd13C12C<-mean(curr.df$d13C12C)
    avgd18O13C<-mean(curr.df$d18O16O/curr.df$d13C12C)
    #avgd18O13C<-mean(curr.df$d18O13C)
    # add to vector
    avgd18O16O.vec<-c(avgd18O16O.vec,avgd18O16O)
    avgd13C12C.vec<-c(avgd13C12C.vec,avgd13C12C)
    avgd18O13C.vec<-c(avgd18O13C.vec,avgd18O13C)
  }

  # make into df
  avgDeltas.df<-as.data.frame(cbind(analysis.vec,avgd18O16O.vec,avgd13C12C.vec,avgd18O13C.vec))
  colnames(avgDeltas.df)<-c("Analysis","avg_d18O16O","avg_d13C12C","avg_d18O13C")
  #head(avgDeltas.df)
  return(avgDeltas.df)
}


# separate full data into list by analysis number
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


# change list into a df
listToDF<-function(dat.list){
  # thread lists into dataframes
  # need total number of elements in the lists
  # sample list
  numSamp<-length(dat.list)
  sampSum<-0
  for(i in seq(1,numSamp)){
    sampIndexLength<-dim(dat.list[[i]])[1]
    sampSum<-sampSum+sampIndexLength
  }
  # initalize df
  numSampCols<-length(colnames(dat.list[[1]]))
  samp.df<-as.data.frame(matrix(rep(NA,numSampCols*sampSum),ncol=numSampCols))
  colnames(samp.df)<-colnames(dat.list[[1]])
  # loop through list and build dfs
  sampIndex<-1
  for(i in seq(1,numSamp)){
    # grab ref vend
    sampDat<-dat.list[[i]]
    numSampPts<-length(sampDat$PeakNr)
    # set indices
    sampInd.vec<-seq(sampIndex,(sampIndex+numSampPts-1))
    # add to df
    samp.df[sampInd.vec,]<-sampDat
    # update ref index
    sampIndex<-sampInd.vec[length(sampInd.vec)]+1
  }
  return(samp.df)
}


# add salt labels
saltsLabs<-function(data.df,salts.vec,relSalts.vec){
  for(i in seq(1,length(salts.vec))){ # loop through salt vector
    ion.vec<-makeSingleIonConcVec(data.df,salts.vec[i])
    data.df[,dim(data.df)[2]+1]<-as.factor(ion.vec)
    colnames(data.df)[length(colnames(data.df))]<-salts.vec[i]
  }

  for(i in seq(1,length(relSalts.vec))){
    ion.vec<-makeIonConcVec(data.df,relSalts.vec[i])
    data.df[,dim(data.df)[2]+1]<-as.factor(ion.vec)
    colnames(data.df)[length(colnames(data.df))]<-relSalts.vec[i]
  }

  for(i in seq(1,length(salts.vec))){ # loop through salt vector
    ion.vec<-make3IonConcVec(data.df,salts.vec[i])
    data.df[,dim(data.df)[2]+1]<-as.factor(ion.vec)
    colnames(data.df)[length(colnames(data.df))]<-paste("rel_",salts.vec[i],sep="")
  }

  return(data.df)
}

# func for saltsLabs
makeSingleIonConcVec<-function(DF,ionConc){
  ionConc.vec<-rep(NA,length(DF$Identifier1))
  ionIndU<-which(grepl(paste(ionConc,"_U",sep=""),DF$Identifier1))
  ionIndL<-which(grepl(paste(ionConc,"_L",sep=""),DF$Identifier1))
  #ionConc.vec[ionIndU]<-paste(ionConc,"_U",sep="")
  #ionConc.vec[ionIndL]<-paste(ionConc,"_L",sep="")
  # replace NAs
  naInd<-which(is.na(ionConc.vec))
  ionConc.vec[naInd]<-paste("No_",ionConc,sep="")
  return(ionConc.vec)
}

# func for saltsLabs
make3IonConcVec<-function(DF,ionConc){
  ionConc.vec<-rep(NA,length(DF$Identifier1))
  ionIndU<-which(grepl(paste(ionConc,"_U",sep=""),DF$Identifier1))
  ionIndL<-which(grepl(paste(ionConc,"_L",sep=""),DF$Identifier1))
  ionConc.vec[ionIndU]<-paste(ionConc,"_U",sep="")
  ionConc.vec[ionIndL]<-paste(ionConc,"_L",sep="")
  # replace NAs
  naInd<-which(is.na(ionConc.vec))
  ionConc.vec[naInd]<-paste("No_",ionConc,sep="")
  return(ionConc.vec)
}

# func for saltsLabs
makeIonConcVec<-function(DF,ionConc){
  ionConc.vec<-rep(NA,length(DF$Identifier1))
  ionInd<-which(grepl(ionConc,DF$Identifier1))
  ionConc.vec[ionInd]<-ionConc
  #replace NA with "No_MgSO4_L"
  naInd<-which(is.na(ionConc.vec))
  ionConc.vec[naInd]<-paste("No_",ionConc,sep="")
  return(ionConc.vec)
}



# salt content
saltContentLab<-function(data.df){
  id1.vec<-data.df$Identifier1
  split<-strsplit(id1.vec,"_")
  saltLab.vec<-c()
  for(i in seq(1,length(split))){
    currSplit<-split[[i]]
    grOne<-c()
    for(j in seq(1,length(currSplit))){
      inSplit<-currSplit[j]
      if(nchar(currSplit[j])>1){
        grOne<-c(grOne,inSplit)
      }
    }
    if(length(grOne)==2){
      saltLab<-paste(grOne[1],"_",grOne[2],sep="")
    }else if(length(grOne)==3){
      saltLab<-paste(grOne[1],"_",grOne[2],"_",grOne[3],sep="")
    }else{
      if(grOne%in%c("H1","L1","LW")){
        grOne<-"no_salt"
      }
      saltLab<-grOne
    }
    saltLab.vec<-c(saltLab.vec,saltLab)
  }
  data.df$saltContent<-saltLab.vec
  return(data.df)
}

## TODO: pH, ionic strength, biotic