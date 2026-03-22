
### Install packages
#install.packages("devtools")
library(devtools)
#install_github("insilico/QCIRMS") 
library(QCIRMS)
library(dplyr)

### Install dependencies
#install.packages(c('dplyr', 'pracma'))
#install_github('isoverse/isoreader')
library(isoreader)
library(icesTAF) # mkdir

# get dxf file names
#dxf_file <- "0_1 m Na2SO4_.dxf" #DXF filename
#file_path = "/Users/vdapoian/FLaRe/QCIRMS/extdata/dxf_files/GSFC_Europa_Experiments/220318_1%CO2_Brines/"
#dataPath <- system.file(file_path, dxf_file, package = "QCIRMS") #DXF file path

SCRIPT_DIR <- "C:/Users/brett-mckinney/Documents/qcirms_tests"
# Relative directory to data
data_dir_name <- "220511_2% CO2 NaCl, Na2SO4, NaHCO3 _Brines"        ## Switch data 
DATA_PATH <- paste(SCRIPT_DIR,"/",data_dir_name,"/",sep="")
PLOTS_PATH <- paste(SCRIPT_DIR,"/QCresults/chromat_plots/",data_dir_name,sep="")
mkdir(PLOTS_PATH)

#### Look for ".dxf"
#setwd(file_path)
#fileNames<-all_dxf_files(path=file_path) # looks in current wd if path not specified 
#length(fileNames)

# plot all raw data using code that writes pdfs in the specified working directory
#generic_plot_all_raw(dxf_file)

dxf_file <- "0_1 m Na2SO4_.dxf" #DXF filename

raw_df <- raw_data(file=dxf_file,path=DATA_PATH)

generic_raw_plot(raw_df,
                 path=PLOTS_PATH, # use path for plots defined before
                 write_pdf=T,
                 title=dxf_file,
                 pdf_name=paste(dxf_file,".pdf",sep=""))
dev.off()


