# Statistics

area_fwhm <- function(matrix3_CopmIDs_all_standard) {
  #returns the area of MZ window (reims) per standard, based on integer (sum over I) for full width at halfmax I
  
  areas <- NULL
  columnnr <- 1
  for (column in matrix3_CopmIDs_all_standard){
    if(columnnr>1){#skip mz col
      
      max_info <- matrix3_CopmIDs_all_standard[which.max(matrix3_CopmIDs_all_standard[,columnnr]),] #to 1 row of df which is max in col x
      max_I <- max_info[columnnr]
      max_mz <- as.numeric(max_info[1])
      halfmaxI <- as.numeric(max_I/2)
      
      df_abovehalfmax <- matrix3_CopmIDs_all_standard[which(matrix3_CopmIDs_all_standard[,columnnr] >= halfmaxI),]
      area <- sum(df_abovehalfmax[,columnnr]) #integreren: sum over area = sumI_width_halfmax
      
      areas <- rbind(areas, area)
    }
    columnnr <- columnnr+1 
  }  
  
  rownames(areas) <- colnames(matrix3_CopmIDs_all_standard[2:length(matrix3_CopmIDs_all_standard)])
  colnames(areas) <- name_standard
  
  return(areas)  #df with area column per sample
}

area_mz <- function(matrix3_CopmIDs_all_standard) {
  #returns the average MZ of peak (reims) per standard
  
  mzs <- NULL
  columnnr <- 1
  for (column in matrix3_CopmIDs_all_standard){
    if(columnnr>1){#skip mz col
      
      max_info <- matrix3_CopmIDs_all_standard[which.max(matrix3_CopmIDs_all_standard[,columnnr]),] #to 1 row of df which is max in col x
      max_mz <- as.numeric(max_info[1])
      
      mzs <- rbind(mzs, max_mz)
    }
    columnnr <- columnnr+1 
  }  
  
  average_mz <- sum(mzs)/length(mzs)
  average_mz <- as.data.frame(average_mz)
  colnames(average_mz) <- name_standard
  
  return(average_mz)  #df with area column per sample
}

linear_eq <- function(x,y) {
  #calculate equation (linear regression) and corresponding R2 value for on plot
  m <- lm(y ~ x -1) #no negative pred conc allowed
  eq <- as.character(
    as.expression(
      substitute(italic(y) == a %.% italic(x) + b*","~~italic(R)^2~"="~r2,
                 list(a = format(round(coef(m)[1], digits=3),nsmall=3), #rico
                      b = 0, #format(round(coef(m)[1], digits=3),nsmall=3), #intercept
                      r2 = format(round(summary(m)$r.squared, digits = 3), nsmall=3)
                      )
      )
    )
  )
  return(eq)
}

Shapiro_Wilk_test <- function(sampleMatrix) {
  #returns the number of normal distribution per variable 
  
  #extra step because intensity in some columns are identical eg. 0
  nrow_matrix <- nrow(sampleMatrix)
  ncol_matrix <- ncol(sampleMatrix)
  random_matrix <- replicate(ncol_matrix, runif(nrow_matrix)*10E-10) 
  sampleMatrix <- sampleMatrix + random_matrix
  
  shapiro20 <- sapply(sampleMatrix, shapiro.test)   #Apply Shapiro Wilk test on each column (=metabolite)
  shapirop20 <- shapiro20[2,]                                     #Retrieve p-values
  shapirop20 <- p.adjust(shapirop20,method="holm",n=length(shapirop20))
  numbernormal20 <-sum(shapirop20 > 0.05)       # chosen alpha value = 0.05
  return(numbernormal20)  #Count number of columns which are normally distributed
}

test_Shapiro_Wilk_test <- function() {
  sampleMatrix <- 'carol-sampleMetadata_short.txt'
  shapiro20 <- sapply(sampleMatrix, shapiro.test)   #Apply Shapiro Wilk test on each column (=metabolite)
  shapirop20 <- shapiro20[2,]                                     #Retrieve p-values
  shapirop20 <- p.adjust(shapirop20,method="holm",n=length(shapirop20))
  numbernormal20 <-sum(shapirop20 > 0.05) 
  
  stopifnot(numbernormal20 == 8 | numbernormal20 == 16) #8 is samples+qc, 16=qc
}

normalize_with_average_QCs <- function(df1, df2){
  #for the intensities of each variable (metabolite): 
  #devides the samples with the mean intenstity of all QC
    #df1 <- samples_metadata
    #df2 <- QC_metadata
  
  df2_matrix <- subset(df2, select= -c(3:COLLUMN_NR_START_VARIABLES-1))
  df2_matrix <- df2_matrix[ ,-1]
  means_QC <- sapply(df2_matrix, mean)  
  
  df1_matrix <- subset(df1, select= -c(3:COLLUMN_NR_START_VARIABLES-1))
  df1_matrix <- df1_matrix[ ,-1]
  
  nrsamples_df1 <- nrow(df1_matrix)
  df_means <- matrix(means_QC,nrow=nrsamples_df1,ncol=length(means_QC),byrow=TRUE)

  normalized_samples_matrix <- df1_matrix/df_means

  # make normalized df
  normalized_samples_metadata <- normalized_samples_matrix
  normalized_samples_metadata$SampleName <- df1$SampleName
  
  #remove inf and NA
  is.na(normalized_samples_metadata) <- sapply(normalized_samples_metadata, is.infinite)
  normalized_samples_metadata[is.na(normalized_samples_metadata)] <- 0

  return(normalized_samples_metadata)
}

normalize_with_IQCs <- function(df){
  #for the intensities of each variable (metabolite): 
  #devides the samples with the mean intenstity of IQC present below the samples
    #df <- samples_IQC_metadata
  
  #sort based on order run samples
  df <- df[order(df$Order),]
  
  #in case of blank/standards in between sampels/QCs, so orders skips number: reset order whitin funtion so no skips
  df$Order <- c(1:nrow(df))
  
  #extra steps for creation of subdfs
  df$Type2 <- df$Type
  df$Type <- as.character(df$Type)
  df$Type2 <- as.character(df$Type2)
  df['Type2'] <- c(NA, head(df$Type2, -1))
  row <- 2 #! first is NA because shift
  df[,"boolean"] <- NA
  collumn_nr_Type2 <-  ncol(df)-1

  for (row in 2:nrow(df)) {
    if (df[row, COLLUMN_NR_TYPE] == df[row, collumn_nr_Type2]) {
      df[row,"boolean"] <- TRUE
    }
    if (df[row, COLLUMN_NR_TYPE] != df[row, collumn_nr_Type2]) {
      df[row,"boolean"] <- FALSE
    }   
  }
  
  #creation of subdfs
  df[,"Nr_subdf"] <- NA
  row <- 2 #! first is NA because shift
  dfnr <- 1

  for (row in 2:nrow(df)) {
    if (df[row,"boolean"] == 'TRUE') {
      df[row,"Nr_subdf"] <- dfnr   
    }
    if (df[row,"boolean"] == 'FALSE') {
      dfnr <- dfnr + 1
      df[row,"Nr_subdf"] <- dfnr  
    }
  }
  subdataframes <- split(df, df$Nr_subdf)
  
  #Note: subdataframes[[1]] == subdataframes$`1`
  #add first row again to df (type2: output NA)
  subdf <- 1
  row1 <- df[1,]
  subdataframes[[subdf]][nrow(subdataframes[[subdf]])+1,] <- row1

  #make sure start with samples as first subdf
  if (subdataframes[[subdf]]$Type[1] == "Sample") {
    subdf <- 1
  }
  if (subdataframes[[subdf]]$Type[1] == "IQC" | subdataframes[[subdf]]$Type[1] == "QC") {
    subdf <- 2
  }
  amount_of_subdfs <- length(subdataframes)
  
  #loop normalize for each subdf with samples
  normalized_samples_matrix <- NULL

  for (subdf in seq(from=subdf, to=amount_of_subdfs-1, by=2)) {    
    df2 <- subset(subdataframes[[subdf+1]], select = -c(Type2, boolean, Nr_subdf))
    df2_matrix <- subset(df2, select= -c(3:COLLUMN_NR_START_VARIABLES-1))
    df2_matrix <- df2_matrix[ ,-1]
    means_IQC <- sapply(df2_matrix, mean)

    df1 <- subset(subdataframes[[subdf]], select = -c(Type2, boolean, Nr_subdf))
    df1_matrix <- subset(df1, select= -c(3:COLLUMN_NR_START_VARIABLES-1))
    df1_matrix <- df1_matrix[ ,-1]

    nrsamples_df1 <- nrow(df1_matrix)
    df_means <- matrix(means_IQC,nrow=nrsamples_df1,ncol=length(means_IQC),byrow=TRUE)

    normalized_samples_df1 <- df1_matrix/df_means
    
    normalized_samples_matrix <- rbind(normalized_samples_matrix, normalized_samples_df1)
  }
  
  #make normalized df
  samples_metadata <- df[df$Type == 'Sample',] 
  samples_metadata <- samples_metadata[order(samples_metadata$Order),]

  #Note: works only because before no rownames defined, default rownames are same as ordernr
  normalized_samples_matrix <- normalized_samples_matrix[order(as.numeric(row.names(normalized_samples_matrix))),]
  normalized_samples_metadata <- normalized_samples_matrix
  normalized_samples_metadata$SampleName <- samples_metadata$SampleName
  
  #remove inf and NA
  is.na(normalized_samples_metadata) <- sapply(normalized_samples_metadata, is.infinite)
  normalized_samples_metadata[is.na(normalized_samples_metadata)] <- 0
   
  return(normalized_samples_metadata)
}

normalize_QC_RLSC <- function(sampleMetadata){
  #Normalize with quality control sample based robust LOESS 
  #(locally estimated scatterplot smoothing) signal correction (QC-RLSC) method 
  #https://cran.r-project.org/web/packages/NormalizeMets/vignettes/NormalizeMets_vignette.html
  library(NormalizeMets)
	
  sampleMetadata <- sampleMetadata[order(sampleMetadata$Order),]
  
  #make featuredata
  samples_IQC_metadata <- sampleMetadata[sampleMetadata$Type == 'QC' | sampleMetadata$Type == 'IQC' | sampleMetadata$Type == 'Sample',] #only samples and QC remain because 1st and last need to be QC
  sampleMetadata <- samples_IQC_metadata #for name below don't need change
  
  #need QC before samples for loess check, st in subset samples/qcs check here
  if(sampleMetadata$Type[1] != "IQC" & sampleMetadata$Type[1] != "QC"){
    stop("Error: LOESS algorithm needs to start with QC/IQC before samples.")
  }

  samplematrix <- from_df_to_matrix(samples_IQC_metadata)
  featuredata <- samplematrix
  rownames(featuredata) <- sampleMetadata$SampleName

  #make sampledata 
  sampledata <- NULL
  sampledata$batch <- sampleMetadata$Batch
  
  for (sample_row in 1:nrow(sampleMetadata)){
  	if (sampleMetadata$Type[sample_row] == "QC" | sampleMetadata$Type[sample_row] == "IQC"){
      sampledata$class[sample_row] <- 0
  	}
  }
  for (sample_row in 1:nrow(sampleMetadata)){
    if (sampleMetadata$Type[sample_row] == "Sample"){
	   sampledata$class[sample_row] <- 1
	 }
  }
  for (sample_row in 1:nrow(sampleMetadata)){
    if (sampleMetadata$Type[sample_row] == "STD" | sampleMetadata$Type[sample_row] == "ISTD"){
	    sampledata$class[sample_row] <- 2
	 }
   }
  for (sample_row in 1:nrow(sampleMetadata)){
    if (sampleMetadata$Type[sample_row] == "Blank"){
	    sampledata$class[sample_row] <- 3
	}
  }
  sampledata$order <- sampleMetadata$Order
  sampledata <- as.data.frame(sampledata)
  rownames(sampledata) <- sampleMetadata$SampleName
	
  #make sure order starts with 1
  sampledata$order <- sapply(sampledata$order, function(x) (x - (ORDER_NR_OF_FIRST_QC-1)))

  #extra step because intensity in some columns are identical eg. 0
  nrow_matrix <- nrow(featuredata)
  ncol_matrix <- ncol(featuredata)
  random_matrix <- replicate(ncol_matrix, runif(nrow_matrix)*10E-10) 
  featuredata <- featuredata + random_matrix
  
  #run LOESS function
  tryCatch(normalized_samples_metadata <- NormQcsamples(featuredata = featuredata, sampledata = sampledata,span = SPAN, lg=FALSE), 
           error=function(cond) print("ERROR: Span is to small, try setting span = 1.") ) #lg = log transformed NOT needed auto so false
  normalized_samples_matrix <- normalized_samples_metadata[["featuredata"]] #save normalised variables
  normalized_samples_matrix <- as.data.frame(normalized_samples_matrix)

  # make df w/o info with SampleName as last collumn
  normalized_samples_matrix$SampleName <- sampleMetadata$SampleName #in rownames feat opgeslagen

  #remove inf and NA (ev verwijder als niet nodig)
  if(nrow(normalized_samples_matrix) != 0){
    is.na(normalized_samples_matrix) <- sapply(normalized_samples_matrix, is.infinite)
    normalized_samples_matrix[is.na(normalized_samples_matrix)] <- 0
  }
  return(normalized_samples_matrix)
}

normalize_not <- function(samples_metadata){
  normalized_samples_metadata <- samples_metadata 
    
  # make df w/o info with SampleName as last collumn
  normalized_samples_matrix <- from_df_to_matrix(normalized_samples_metadata) 
  normalized_samples_matrix$SampleName <- normalized_samples_metadata$SampleName

  #remove inf and NA
  if(nrow(normalized_samples_matrix) != 0){
      is.na(normalized_samples_matrix) <- sapply(normalized_samples_matrix, is.infinite)
      normalized_samples_matrix[is.na(normalized_samples_matrix)] <- 0
  }  
  return(normalized_samples_matrix)
}
  
log_transformation <- function(normalizedMetadata){
  # Log transform: (1+x, base=exp(1))
  #natural log of value+1 so works accurately for negative vanues -1<x<0
  lognormalizedMetadata <- normalizedMetadata
  lognormalizedMetadata[,COLLUMN_NR_START_VARIABLES:length(normalizedMetadata)] <- sapply(normalizedMetadata[,COLLUMN_NR_START_VARIABLES:length(normalizedMetadata)], log1p)   #Transform columns
  
  return(lognormalizedMetadata)
}

Pareto_scaling <- function(lognormalizedMetadata){
  # Scaling: pareto: (observation - column mean)/sqrt(stdev column)
  #synonym standard scaler maar met de wortel
  x.centered <- lognormalizedMetadata
  x.centered[,COLLUMN_NR_START_VARIABLES:length(lognormalizedMetadata)] <- apply(lognormalizedMetadata[,COLLUMN_NR_START_VARIABLES:length(lognormalizedMetadata)], 2, function(x) x - mean(x))
  scalednormalizedMetadata <- x.centered
  scalednormalizedMetadata[,COLLUMN_NR_START_VARIABLES:length(lognormalizedMetadata)] <- apply(x.centered[,COLLUMN_NR_START_VARIABLES:length(lognormalizedMetadata)], 2, function(x) x/sqrt(sd(x)))
  
  return(scalednormalizedMetadata)
} 

critical.r <- function( n, alpha = .05 ) {
  df <- n - 2
  critical.t <- qt(alpha/2, df, lower.tail = F)
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return(critical.r)
}

correlation_matrix_samples <- function(scaledmatrix_samples1, scaledmatrix_samples2) {
  ## Correlation coeficient calcualtion
  
  #optional 2nd arg
  if(missing(scaledmatrix_samples2)) {
    correlation_scaledmatrix <- round(cor(scaledmatrix_samples1, method = "spearman"),3)
  } else {
    correlation_scaledmatrix <- round(cor(scaledmatrix_samples1, scaledmatrix_samples2, method = "spearman"),3)
  }
  
  return(correlation_scaledmatrix)
}  

Pvalues_correlation_matrixFDR <- function(scaledmatrix_samples1, scaledmatrix_samples2) {
  ## Correlation coeficient adjusted p-values calculation
  suppressMessages(library(psych))
  
  #optional 2nd arg
  if(missing(scaledmatrix_samples2)) {
    adj_correlation_scaledmatrix <- corr.test(scaledmatrix_samples1, use="complete", method="spearman", adjust="fdr", alpha=.05)
    #matrix$r same as using corr in default func but slower
  } else {
    adj_correlation_scaledmatrix <- corr.test(scaledmatrix_samples1, scaledmatrix_samples2, use="complete", method="spearman", adjust="fdr", alpha=.05)
    #matrix$r same as using corr in default func but slower
  }
  
  Pvalues_correlation_scaledmatrix <- round(adj_correlation_scaledmatrix$p, 3) #extract FDR adjusted p-values
  
  detach(package:psych) #issue mixomics, ggplot possibly
  return(Pvalues_correlation_scaledmatrix)
}

correlation_analysis_loop <- function(scaledmatrix_samples1, scaledmatrix_samples2){  #merged_df
  #calculates correlation coeffic (Spearman, FDR) per row iteration
  
  suppressMessages(library(psych))
  
  #only on merged, so both intr and inter calc, need 2dfs for func
  A.data <- scaledmatrix_samples1 #(merged_df)#[1:10, 1:10])
  B.data <- scaledmatrix_samples2 #(merged_df)#[1:10, 1:10])
  
  # create a presized data frame
  Arows <- nrow(A.data)
  Brows <- nrow(B.data)
  correlation.result <- data.frame(t(matrix(c(0.0, 0.0, "foo", "foo"), nrow=4, ncol=Arows * Brows)), stringsAsFactors=FALSE)
  colnames(correlation.result) <- c("rho", "adj_p_value", "CompID1", "CompID2")
  
  # assign values to that presized data frame 
  Acols <- ncol(A.data)
  Bcols <- ncol(B.data)
  for (i in 1:Arows) {
    # get A data
    Ai <- as.numeric(A.data[i, 1:Acols])
    Aname <- rownames(A.data)[i]
    for (j in 1:Brows) {
      # calculate row offset in presized data frame
      ij <- (i - 1) * Brows + j
      # get B data
      Bj <- as.numeric(B.data[j, 1:Bcols])
      Bname <- rownames(B.data)[j]
      # run test
      correln <- cor.test(Ai, Bj, use="complete", method="spearman", adjust="fdr", alpha=.05)
      rho <- round(correln$estimate,3)
      p.val <- round(correln$p.value,3)
      source <- paste(Aname)
      metab <- paste(Bname)
      # assign value to row of result data frame
      correlation.result[ij,] <- c(rho, p.val, source, metab)
    }
  }
  detach(package:psych) #issue mixomics, ggplot possibly
  
  return (correlation.result)
}


