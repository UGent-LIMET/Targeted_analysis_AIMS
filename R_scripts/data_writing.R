# data writing

append_result_to_report <- function(result, outputFile){
  # write line results to report.txt file, with append
  write(result, file=outputFile, append=TRUE)
 
  stopifnot(length(result) > 0)
  
  return()
}

append_dataframe_to_report <- function(result, outputFile){
  # write multiplelines as df to report.txt file, with append
  write.table(result, file=outputFile, append = TRUE, sep='\t', col.names = FALSE, row.names = FALSE)
  stopifnot(length(result) > 0)
  return()
}

write_matrix_as_txt_file <- function(dataframe, outputFile){
  # write df (w/o the 'fake' title header) and w col as title (col1) results to txt file
  write.table(dataframe, file=outputFile, sep ="\t", row.names = TRUE, col.names = FALSE)
 
  stopifnot(nrow(dataframe) >= 1)
  
  return()
}

write_matrixTT_as_txt_file <- function(dataframe, outputFile){
  # write df (with title header)(row1) and col1 results to txt file, but move 1st col so empty top-left cell
  
  # move 1st col so empty top-left cell
  a <- rep("", nrow(dataframe))
  dataframe2 <- cbind(dataframe, a)
  colnames(dataframe2) <- c("", colnames(dataframe))
  write.table(dataframe2, file=outputFile, sep ="\t", row.names = TRUE, col.names = TRUE)
  
  stopifnot(nrow(dataframe2) >= 1)
  
  return()
}

write_dataframe_as_txt_file <- function(dataframe, outputFile){
  # write df (with title header)(row1) results to txt file
  write.table(dataframe, file=outputFile, sep ="\t", row.names = FALSE, col.names = TRUE)
 
  stopifnot(nrow(dataframe) >= 1)
  
  return()
}

write_dataframeTT_as_txt_file <- function(dataframe, outputFile){
  # write df (with title header)(row1) and col1 results to txt file
  write.table(dataframe, file=outputFile, sep ="\t", row.names = TRUE, col.names = TRUE)
 
  stopifnot(nrow(dataframe) >= 1)
  
  return()
}

