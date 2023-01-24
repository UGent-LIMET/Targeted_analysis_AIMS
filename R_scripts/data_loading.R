# All data loading from external sources

load_sampleMetadata <- function(INPUT_SAMPLES){
  # load samplelist from .txt file
  try(sampleMetadata <- read.table(INPUT_SAMPLES, header=TRUE, sep="\t"))
  if(exists("sampleMetadata") == FALSE){
    try(sampleMetadata <- read.table(file=INPUT_SAMPLES, header=TRUE, sep="\t", fileEncoding="UTF-16")) 
  }
  if(exists("sampleMetadata") == FALSE){
    stop('Can not read SM.txt files, check if UTF-16/UTF-8/ASCII format')
  }
  
  stopifnot("SampleName" %in% colnames(sampleMetadata))
  stopifnot(colnames(sampleMetadata)[1] == "SampleName")
  
  
  indx <- as.numeric(substring(sampleMetadata$SampleName[1], 1, 1)) 
  stopifnot(NA %in% indx == FALSE) #variables must start with number
  
  stopifnot("Type" %in% colnames(sampleMetadata))
  stopifnot("Order" %in% colnames(sampleMetadata))
  stopifnot("Batch" %in% colnames(sampleMetadata))
  
  stopifnot(ncol(sampleMetadata) == (COLLUMN_NR_PROJECTION1 + AMOUNT_OF_PROJECTIONS-1)) #projection is last collumn
  
  min_samples <- 2
  stopifnot(nrow(sampleMetadata) >= min_samples) #min amount of samples needed
  
  return(sampleMetadata)
}

test_load_sampleMetadata <- function(){
  fn2 <- 'carol-sampleMetadata_short.txt'
  sampleMetadata <- read.table(fn2, header=TRUE, sep="\t")
  
  stopifnot(nrow(sampleMetadata) == 40)
}

load_variableMetadata <- function(INPUT_VARIABLES){
  # load samplelist from .txt file
  #todo!!!! need distribution min/max to work!!! add tests in read VM file
  variableMetadata <- read.table(INPUT_VARIABLES, header=TRUE, fill=T, sep="\t")
 
  stopifnot(class(variableMetadata$CompID)=="integer") #CompID are integers (eg. 1, 2, ..., n)
  indx <- apply(variableMetadata[,COLLUMN_NR_START_SAMPLES:ncol(variableMetadata)], 2, function(x) any(is.na(x) | is.infinite(x) | is.numeric(x)))
  stopifnot(TRUE %in% indx == TRUE) #intensities variables are float numbers with decimal (eg. 1000000.00), no NA/inf/"1,000,000" values
  stopifnot(min(variableMetadata[,COLLUMN_NR_START_SAMPLES:ncol(variableMetadata)]) >= 0) #no negative or empty numbers allowed
  
  stopifnot("CompID" %in% colnames(variableMetadata))
  stopifnot("MZ" %in% colnames(variableMetadata))
  stopifnot("Time" %in% colnames(variableMetadata))
  stopifnot("isotopes" %in% colnames(variableMetadata))
  stopifnot("adduct" %in% colnames(variableMetadata))
  stopifnot("pcgroup" %in% colnames(variableMetadata))
  stopifnot("name" %in% colnames(variableMetadata))
  stopifnot("fold" %in% colnames(variableMetadata))
  stopifnot("tstat" %in% colnames(variableMetadata))
  stopifnot("pvalue" %in% colnames(variableMetadata))
  stopifnot("mzmed" %in% colnames(variableMetadata))
  stopifnot("mzmin" %in% colnames(variableMetadata))
  stopifnot("mzmax" %in% colnames(variableMetadata))
  stopifnot("rtmed" %in% colnames(variableMetadata))
  stopifnot("rtmin" %in% colnames(variableMetadata))
  stopifnot("rtmax" %in% colnames(variableMetadata))
  stopifnot("npeaks" %in% colnames(variableMetadata))
  stopifnot("bio" %in% colnames(variableMetadata))
  stopifnot("blank" %in% colnames(variableMetadata))
  
  min_samples <- 2
  stopifnot(nrow(variableMetadata) >= min_samples) #more then 10 samples in analysis
  
  return(variableMetadata)
}

load_targetedMetadata <- function(INPUT_STANDARDS){
  # load standards list from .txt file
  targetedMetadata <- read.table(INPUT_STANDARDS, header=TRUE, sep="\t")
  
  #all columns present
  stopifnot("STD_ID" %in% colnames(targetedMetadata)) 
  stopifnot("MetaboliteName" %in% colnames(targetedMetadata)) 
  stopifnot("IonisationMode" %in% colnames(targetedMetadata)) 
  stopifnot("MZ" %in% colnames(targetedMetadata))
  stopifnot("Time" %in% colnames(targetedMetadata)) #need 1 time (can have others in next cols as well and this col = average)

  #MZ numeric 
  indx <- any(is.na(targetedMetadata$MZ) | is.infinite(targetedMetadata$MZ) | is.numeric(targetedMetadata$MZ))
  stopifnot(TRUE %in% indx == TRUE) #intensities variables are float numbers with decimal (eg. 1000000.00), no NA/inf/"1,000,000" values
  stopifnot(min(targetedMetadata$MZ) >= 0) #no negative numbers allowed
  #stop if there are rows not defines as "+/-" in ionisation mode  (eg. ?, +Na, NA)
  stopifnot(nrow(targetedMetadata) == nrow(targetedMetadata[targetedMetadata$IonisationMode == '+' | targetedMetadata$IonisationMode == '-',]))
  
  #time num.int
  #stopifnot(class(targetedMetadata$Time)=="numeric" | class(targetedMetadata$Time)=="integer") #if empty for reims: int OR logi..., else numbr
  
  min_samples <- 1
  stopifnot(nrow(targetedMetadata) > min_samples) #min 2stds really needed for not wrong merge...
  
  return(targetedMetadata)
}

load_matrix3_CopmIDs_all <- function(){
  #loads load_matrix3_CopmIDs_all (present in scans folder after pre-processing REIMS)
  
  namedf <- paste0(path_data_out, "/scans/matrix3_CopmIDs_all.txt") 
  if(file.exists(namedf) == FALSE){
    stop('Can not perform visualisation, pre-processing needs to be performed in R pipeline')
  }
  
  matrix3_CopmIDs_all <- read.table(file=namedf, header=TRUE, sep="\t")
  
  return(matrix3_CopmIDs_all)
}

load_large_txt_file <- function(file_){
  #loads large txt file (eg. converted from thermo or waters raw files)
  
  library(data.table) #fread laden
  
  try(scans <- fread(file=file_, header=TRUE, sep="\n"))  #enter seps, header= number of scans
  if(exists("scans") == FALSE){
    try(scans <- fread(file=file_, header=TRUE, sep="\n", fileEncoding="UTF-16")) 
  }
  if(exists("scans") == FALSE){
    stop('Can not read raw txt file, check if UTF-16/UTF-8/ASCII format')
  }
  
  return(scans)
}

load_large_df <- function(file_){
  #loads large df (eg. converted from termo txt file)
  
  library(data.table) #fread laden
  
  try(df_file <- fread(file=file_, header=TRUE, sep="\t"))  #enter seps, header= number of scans
  if(exists("df_file") == FALSE){
    try(df_file <- fread(file=file_, header=TRUE, sep="\t", fileEncoding="UTF-16")) 
  }
  if(exists("df_file") == FALSE){
    stop('Can not read raw txt file, check if UTF-16/UTF-8/ASCII format')
  }
  
  return(df_file)
}

load_correlationMetadata <- function(INPUT_VARIABLES1){
  # load samplelist from .txt file
  try(correlationMetadata <- read.table(INPUT_VARIABLES1, header=TRUE, sep="\t"))
  if(exists("correlationMetadata") == FALSE){
    try(correlationMetadata <- read.table(file=INPUT_VARIABLES1, header=TRUE, sep="\t", fileEncoding="latin1")) #latin1 on ibof, ok on r reims
  }
  if(exists("correlationMetadata") == FALSE){
    stop('Can not read CM.txt files, check if UTF-16/UTF-8/ASCII format')
  }
  
  stopifnot("SampleName" %in% colnames(correlationMetadata))
  stopifnot(colnames(correlationMetadata)[1] == "SampleName")
  
  min_samples <- 2
  stopifnot(nrow(correlationMetadata) >= min_samples) #min amount of samples needed
  
  return(correlationMetadata)
}

load_variableMetadata_other <- function(INPUT_VARIABLES, COLLUMN_NR_START_SAMPLES){
  #less strict than above, also other variables than metab allowed
  #also neg numbers and NAs allowed
  
  variableMetadata <- read.table(INPUT_VARIABLES, header=TRUE, fill=T, sep="\t")
  
  #compid needs to be present but not int
  #stopifnot(class(variableMetadata$CompID)=="integer") #CompID are integers (eg. 1, 2, ..., n)
  
  indx <- apply(variableMetadata[,COLLUMN_NR_START_SAMPLES:ncol(variableMetadata)], 2, function(x) any(is.na(x) | is.infinite(x) | is.numeric(x)))
  
  #now negative numbers and NAs ARE allowed:
  #stopifnot(TRUE %in% indx == TRUE) #intensities variables are float numbers with decimal (eg. 1000000.00), no NA/inf/"1,000,000" values
  #stopifnot(min(variableMetadata[,COLLUMN_NR_START_SAMPLES:ncol(variableMetadata)]) >= 0) #no negative or empty numbers allowed
  
  stopifnot("CompID" %in% colnames(variableMetadata))

  min_samples <- 2
  stopifnot(nrow(variableMetadata) >= min_samples) #more then 10 samples in analysis
  
  return(variableMetadata)
}
