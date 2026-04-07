library(dplyr)
library(isoreader)
library(QCIRMS)

####### new data merging
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ts_path <- "./feature_extraction/ts_feature_extraction_results/"
ms_path <- "./feature_extraction/ms_feature_extraction_results/"


#####

# read data in as lists - each elem is a file
ts_files <- list.files(path=ts_path)
head(ts_files)
length(ts_files)
ms_files <- list.files(path=ms_path)
head(ms_files)
length(ms_files)

ts.list <- list()
for(i in seq(1,length(ts_files))){
  curr_dat <- read.csv(paste(ts_path, ts_files[i], sep=""))
  ts.list[[i]] <- curr_dat
}
length(ts.list)
head(ts.list[[1]])
dim(ts.list[[26]])
# [1] 72 54

ts.df <- do.call("rbind", ts.list)
dim(ts.df)
# [1] 1308   54

ms.list <- list()
for(i in seq(1,length(ms_files))){
  curr_dat <- read.csv(paste(ms_path, ms_files[i], sep=""))
  ms.list[[i]] <- curr_dat
}
length(ms.list)
head(ms.list[[1]])
dim(ms.list[[26]])
# [1] 72 67

ms.df <- do.call("rbind", ms.list)
dim(ms.df)
# [1] 1309   67


## find which experiment is missing in the ts list
ms_files.26 <- ms.list[[26]]$fileId
ms_an.26 <- ms.list[[26]]$Analysis
length(ms_an.26)
# [1] 35
length(unique(ms_an.26))
# [1] 34
ms_an.26[which(duplicated(ms_an.26))]
# [1] 15897

ms_files.26[which(ms_an.26 == "15897")]
# [1] "L1_(10).dxf" "L1_(9).dxf"
# data looks identical, so just remove the one of the duplicates

ms.list[[26]][which(ms_an.26 == "15897"),]
ms_files[26]

ts_an.26 <- ts.list[[26]]$Analysis




#### code from previous processing
# remove replicate ms experiments
# remove standards experiments
msfeat_norep <- ms.df %>%
                group_by(Analysis)     %>% # for each Analysis group
                slice(n())             %>% # grab last row of unique Analysis
                # rm rows that are standards:
                #filter(!(Identifier1 %in% c("H1", "L1", "LW") )) %>%
                filter(!(grepl('^H1', Identifier1))) %>%
                filter(!(grepl('^L1', Identifier1))) %>%
                filter(!(grepl('^LW', Identifier1))) %>%
                filter(!grepl("CO2", Identifier1)) %>% 
                data.frame()
dim(msfeat_norep)  # [1] 847  67, unique Analysis numbers and removed standards


######## Merge
# merge does the intersection
ms_ts_merge_all_cols <- merge(msfeat_norep, ts.df, by="Analysis")

######## Filter and create the mass-spec + time series data
# columns to keep in merged data
# keep ts features, includes "Analysis" -- probably keep Analysis for later merging
# keep ms features that start with "avg"
ms_ts_data <- ms_ts_merge_all_cols %>%
  select(colnames(ts.df), starts_with("avg")) # add biotic to data?
# remove 0 variation data
sd_filter <- apply(ms_ts_data,2,function(x) {sd(x)})
ms_ts_data <- ms_ts_data[, sd_filter!=0]

dim(ms_ts_data)  # [1] 853  82

######## Create the meta data
meta_data_features <- c("Analysis", "fileId", "Identifier1", "Preparation","DateTime")
meta_data <- ms_ts_merge_all_cols %>%
  select(one_of(meta_data_features))

dim(meta_data)
colnames(meta_data)
colnames(ms_ts_data)

write.csv(ms_ts_data,file="./feature_extraction/new_ms_ts_data.csv", row.names=F)
write.csv(meta_data,file="./feature_extraction/new_meta_data.csv", row.names=F)

###############
### add ML labels
## TODO: pH and ionic strength?
