# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part I bis: targeted analysis REIMS


##########R Pipeline - Part I bis: targeted analysis##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part I bis: targeted analysis REIMS - start!")
# Part I bis: targeted analysis REIMS


## data_loading
setwd(path_data_in)
if (INPUT_VARIABLES == VARIABLEMETADATA_FROM_PIPELINE){
  INPUT_VARIABLES <- paste(name_project, '_variableMetadata.txt', sep="")
  if(exists("COLLUMN_NR_START_SAMPLES") == FALSE){ 		#if after merge, will be given value 21, so do not change
    COLLUMN_NR_START_SAMPLES <- 20  #always the same if auto == 20
  }
}

variableMetadata <- load_variableMetadata(INPUT_VARIABLES) 
#note: automatically adds X to samples at this point

targetedMetadata <- load_targetedMetadata(INPUT_STANDARDS) 
#always same since in read_me

#specific for reims TM:
indx <- min(targetedMetadata$MZ) < 50 | max(targetedMetadata$MZ) > 1200
if(TRUE %in% indx == TRUE){
  stop("remove MZs with values < 50 and > 1200 Da from TargetedMetadata")
}


if(VISUALISATION_MZ_WINDOW_PPM == VISUALISATION |AREA == CALCULATE_AREA){
  matrix3_CopmIDs_all <- load_matrix3_CopmIDs_all()
}


## set directory to output
setwd(path_data_out)


## retain only correct ionis mode + use coorect mz values based on ionisation mode
if (POLARITY == "positive") {
  retained_targetedMetadata <- targetedMetadata[targetedMetadata$IonisationMode == "+", ]
}
if (POLARITY == "negative") {
  retained_targetedMetadata <- targetedMetadata[targetedMetadata$IonisationMode == "-", ]
}


##for each variable, search als collect in output pre-processing
variableMetadata_standards <- NULL
matrix3_CopmIDs_all_standards <- NULL

number_of_samples <- ncol(variableMetadata) - COLLUMN_NR_START_SAMPLES+1
areas_df <- data.frame(matrix(vector(), nrow=number_of_samples, ncol=(1), #entry for each sample
                              dimnames=list(c(), "SampleName")),
                       stringsAsFactors=F)
MZs <- data.frame(matrix(vector(), nrow=1, ncol=(1), #entry: 1row with average mz
                              dimnames=list(c(), "MZ")),
                       stringsAsFactors=F)

ppm <- PPM #3 for Q-exactutive (lipids), 5 for exacutive and Q-exactutive (polar), 50 for REIMS

for(row in 1:nrow(retained_targetedMetadata)){  
  
  #row <- 13
  
  standard <- as.numeric(as.character(retained_targetedMetadata$MZ[row])) #stomme truk nodig voor factors..
  name_standard <- as.character(retained_targetedMetadata$MetaboliteName[row]) #Add name standard to variableMetadata
  standard_down <- standard - ((ppm*standard)/1000000) 
  standard_up <- standard + ((ppm*standard)/1000000)
  variableMetadata_standard <- variableMetadata[variableMetadata$MZ >= standard_down & variableMetadata$MZ <= standard_up, ]

  
  ##PLOT MZ dimension
  if(VISUALISATION_MZ_WINDOW_PPM == VISUALISATION){
    matrix3_CopmIDs_all_standard <- matrix3_CopmIDs_all[matrix3_CopmIDs_all$MZ >= standard_down & matrix3_CopmIDs_all$MZ <= standard_up, ]
    #make plot for manual evaluation
    #plot(matrix3_CopmIDs_all_standard[,c(1,2)],type="b") #mz vs I of 1s sample
    
    #also name to write boxplot.png (so no "/" eg)
    library("stringr") 
    name_standardWOsymbol <- str_replace_all(name_standard, "[^[:alnum:]]", "")    # Delete non-alphanumeric
    
    #all samples w ggplot
    library("reshape")
    data_melt <- melt(matrix3_CopmIDs_all_standard, id = c("MZ")) 
    
    if(number_of_samples <= 15){
      #plot wo legend
      MZdimension_plot <- plot_MZplane(data_melt, ppm)
    
      png(paste(name_project, "_MZdimension_std_", name_standardWOsymbol, ".png", sep=""), width=10, height=5, units="in", res=150)
      plot(MZdimension_plot)
      dev.off()
    }
    if(number_of_samples > 15){
     #plot wo legend
      MZdimension_plot <- plot_MZplaneWOlegend(data_melt, ppm)
      
      png(paste(name_project, "_MZdimension_std_", name_standardWOsymbol, ".png", sep=""), width=7, height=5, units="in", res=150)
      plot(MZdimension_plot)
      dev.off()
      
      #write legend for max 90? samples, also dep length samplenames
      if(row == 1){ #width png 10 = 15*3= 45samples
        MZdimension_plot <- plot_MZplane(data_melt, ppm)
        
        png(paste(name_project, "_MZdimension_LEGEND.png", sep=""), width=20, height=5, units="in", res=150)
        plot(MZdimension_plot)
        dev.off()
      } 
    }
  }
  
  ##calcualte areas
  if(AREA == CALCULATE_AREA){
    matrix3_CopmIDs_all_standard <- matrix3_CopmIDs_all[matrix3_CopmIDs_all$MZ >= standard_down & matrix3_CopmIDs_all$MZ <= standard_up, ]
    
    #area calc (fwhm sum I)
    area_df <- area_fwhm(matrix3_CopmIDs_all_standard)
    area_df <- as.data.frame(area_df)
    
    #report 1 mz, average of mzs per standard
    MZ <- area_mz(matrix3_CopmIDs_all_standard)
    
    areas_df <- cbind(areas_df, area_df)
    MZs <- cbind(MZs, MZ)
  }
  
  ##make df results
  #if no found, continue
  #if >=1 found, write all. ev todo >1: 1select closest to mz
  if(nrow(variableMetadata_standard) >= 1){  
    variableMetadata_standard$MetaboliteName <- name_standard
    variableMetadata_standards <- rbind(variableMetadata_standards, variableMetadata_standard)
  }
}


## write VM small list with only variables found
#add names
variableMetadata_output <- cbind(variableMetadata_standards$CompID, variableMetadata_standards$MetaboliteName) #compIDs subset from full VM after preprocessing
variableMetadata_output <- cbind(variableMetadata_output, variableMetadata_standards[,-1]) #wo compID double
variableMetadata_output <- variableMetadata_output[,-ncol(variableMetadata_output)] #wo info MetaboliteName double
colnames(variableMetadata_output)[1] <- "CompID" #name CompID not correct otherwise
colnames(variableMetadata_output)[2] <- "MetaboliteName" 

#remove "X" before samplename (in P2 load added again but need to be same as VM xcms output for Part2...)
samplenames <- colnames(variableMetadata_output)[COLLUMN_NR_START_SAMPLES+1:number_of_samples]
samplenames <- substr(samplenames, 2,nchar(samplenames))
colnames(variableMetadata_output)[COLLUMN_NR_START_SAMPLES+1:number_of_samples] <- samplenames


setwd(path_data_in) #todo rm comment, to replace input to connect all modules
#write_dataframe_as_txt_file(variableMetadata_output, paste(name_project, '_variableMetadata.txt', sep="")) #input for R_pipeline part II: statistical analysis OR/AND univariate analysis
setwd(path_data_out)
write_dataframe_as_txt_file(variableMetadata_output, paste(name_project, '_variableMetadata_output_targeted_analysis.txt', sep="")) #copy for user


## write VM small list with with areas from M3-based plot
if(AREA == CALCULATE_AREA){
  MZs <- MZs[,2:ncol(MZs)] #rm addit 1st col
  areas_df <- areas_df[,2:ncol(areas_df)] #rm addit 1st col
  standardnames <- colnames(areas_df)
  samplenames <- rownames(areas_df) #with X
  areas_df <- as.data.frame(t(areas_df)) 
  
  number_compIDs <- nrow(areas_df)
  
  variableMetadata_output_bis <- NULL
  variableMetadata_output_bis$CompID <- as.integer(seq(1, number_compIDs, by=1)) #compid name 1,2,.. since cannot subset from original VM preprocessing (from M3 here) <=> STD_ID eg ID001 from lca
  variableMetadata_output_bis$MetaboliteName <- standardnames
  variableMetadata_output_bis$MZ <- as.numeric(MZs) #order never changed so match, ok
  variableMetadata_output_bis$Time <- rep("1",number_compIDs) #put RT in mins
  
  rest_colnames <- c("isotopes",	"adduct",	"pcgroup",	"name",	"fold",	"tstat",	"pvalue",	"mzmed",	"mzmin",	"mzmax",	"rtmed",	"rtmin",	"rtmax",	"npeaks",	"bio",	"blank"	)
  rest_VM <- data.frame(matrix(vector(), nrow=number_compIDs, ncol=16, 
                               dimnames=list(c(), rest_colnames)),
                        stringsAsFactors=F) #empty values for isotopes etc.
  variableMetadata_output_bis <- cbind(variableMetadata_output_bis, rest_VM)
  
  #add sample intensities
  variableMetadata_output_bis <- cbind(variableMetadata_output_bis, areas_df)
  
  #remove "X" before samplename (in P2 load added again but need to be same as VM xcms output for Part2...)
  samplenames <- substr(samplenames, 2,nchar(samplenames))
  colnames(variableMetadata_output_bis) <- c("CompID", "MetaboliteName", "MZ", "Time", rest_colnames, samplenames)
  
  setwd(path_data_in) #todo rm comment, to replace input to connect all modules
  #write_dataframe_as_txt_file(variableMetadata_output, paste(name_project, '_variableMetadata.txt', sep="")) #input for R_pipeline part II: statistical analysis OR/AND univariate analysis
  setwd(path_data_out)
  write_dataframe_as_txt_file(variableMetadata_output_bis, paste(name_project, '_variableMetadata_output_targeted_analysis_bis.txt', sep="")) #copy for user
}

#to work in pipeline, extra col VM so set here
COLLUMN_NR_START_SAMPLES <- COLLUMN_NR_START_SAMPLES + 1  


print("R pipeline - Part I bis: targeted analysis REIMS - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################
