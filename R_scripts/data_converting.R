# data handling

merge_accoding_to_SampleName <- function(df1, df2) {
  # sort both dfs based on SampleName + merge
  df1 <- df1[order(df1$SampleName),]
  df2 <- df2[order(df2$SampleName),]

  stopifnot(nrow(df1) == nrow(df2))
  mergeddf <- merge(df1, df2, by="SampleName", all = T, sort=F)

return(mergeddf)
}    

from_df_to_matrix <- function(df) {
  
  #  # output is numeric df without info (chr/factors/...)
  #  matrix <- subset(df, select= -c(2:(COLLUMN_NR_START_VARIABLES-1)))
  #  matrix <- matrix[ ,-1]
  
  matrix <- df[,COLLUMN_NR_START_VARIABLES:ncol(df)]
  
  return(matrix)
}

merge_two_variableMetadata <- function(df1, df2){
  #df1 <- variableMetadata1 #for testing
  #df2 <- variableMetadata2
  
  # merge two variablemetadata via rbind (eg +/- modi or lipid/polar; or both by calling function twice)
  stopifnot(ncol(df1) == ncol(df2))
  
  mergeddf <- rbind(df1, df2)
  colnames(mergeddf)[colnames(mergeddf) == "CompID"] <- "CompID_merged" #todo check op twice run so no unique colname?
  mergeddf$CompID <- 1:nrow(mergeddf)
  
  #compID as first column
  CompID <- mergeddf$CompID
  all_but_last_cols <- mergeddf[,1:(ncol(mergeddf)-1)]
  mergeddf <- cbind(CompID, all_but_last_cols)
  return(mergeddf)
}

merge_two_variableMetadata_add_samples <- function(df1, df2){
  #df1 <- variableMetadata1 #for testing
  #df2 <- variableMetadata2
  
  #chech not content but already if same amount compIDs
  stopifnot(nrow(df1) == nrow(df2))
  
  #remove redundant info of compIDs in 2nd df
  df2 <- df2[,COLLUMN_NR_START_SAMPLES:ncol(df2)]
  
  # merge two variablemetadata via cbind (eg add new samples in timecource if in vm exact same CompIDs+ordered)
  mergeddf <- cbind(df1, df2)

  return(mergeddf)
}

rtBurned <- function(msFile){
  #preprocessing reims: search where burns
  write(paste(fileNames(msFile), '\\n'), '~/.error.log', append=TRUE)
  threshold <- mean(tic(msFile))
  rtBurns <- which(tic(msFile) >= threshold)

  scansBurned <- split(rtBurns, cumsum(c(TRUE, diff(rtBurns) > 1))) %>% 
    lapply(function(x) x[which.max(tic(msFile)[x])])

  # for each scan search the most intense mz
  mzsBurned <- lapply(scansBurned, function(scan){
    mzs <- mz(msFile[[scan]])
    mzs[which.max(intensity(msFile[[scan]]))]
  })

  lapply(1:length(mzsBurned), function(i)
    rtBurnedByZone(msFile, mzsBurned[[i]], scansBurned[[i]]))
}
  
rtBurnedByZone <- function(msFile, mzBurned, scanBurned){
  #preprocessing reims: search inside burns
  #+/-1 voor exclude aan begin/einde chromatogram 'slechte' stuk
  left <- if(scanBurned > 1) scanBurned -1 else 1
  foundMzBurned <- TRUE
  while(left != 1 & foundMzBurned){
    data <- data.frame(mz = mz(msFile[[left]]), int = intensity(msFile[[left]]))
    foundMzBurned <- mzBurned %in% (data %>% top_n(5, int) %>% pull(mz))
    left <- left - 1
  }
  right <- if(scanBurned < length(msFile)) scanBurned + 1 else length(msFile)
  foundMzBurned <- TRUE
  while(right != length(msFile) & foundMzBurned){
    data <- data.frame(mz = mz(msFile[[right]]), int = intensity(msFile[[right]]))
    foundMzBurned <- mzBurned %in% (data %>% top_n(5, int) %>% pull(mz))
    right <- right + 1
  }
  list(rtmin = rtime(msFile[left]), rtmax = rtime(msFile[right]))
}

make_mz_step <- function(df){
  #function make mz_step correct count until new scan starts
  
  #reuse function iqc, so rename column
  df$Type <- df$scannumber 
  
  #copy column 1 cel lower
  df$Type2 <- df$Type
  df$Type <- as.character(df$Type)
  df$Type2 <- as.character(df$Type2)
  df['Type2'] <- c(NA, head(df$Type2, -1))
  row <- 2 #! first is NA because shift
  
  #add nr into new column
  df$step_mz[1] <- 1
  for (row in 2:nrow(df)) {
    if (df[row, 'Type'] == df[row, 'Type2']) {
      df$step_mz[row] <-  df$step_mz[row-1] + 1
    }
    if (df[row, 'Type'] != df[row, 'Type2']) {
      df$step_mz[row] <- 1 
    }   
  }
  
  #del new columns
  df$Type <- NULL
  df$Type2 <- NULL
  
  return(df)
}
