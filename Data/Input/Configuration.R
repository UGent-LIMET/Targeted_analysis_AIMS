## @knitr INFO
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Configuration


##########Global_settings##########

## @knitr settings
## options
RUN_CODE <- 'run this part of the pipeline'
DONT_RUN_CODE <- 'skip this part of the pipeline, keeps order: pre-processing - targeted analysis - statistical analysis - annotation'

## Adjustments
#Project name:
EXPERIMENT <- 'Targeted_analysis_AIMS' #structured and short, see READ_ME
POLARITY <- "both" #{"positive", "negative"} #needed this format for annotation to work
# file_conversion and pre-processing need to be performed seperate
# only choose "both" when no pre-processing needed (eg. merge both ionisationmodes)
USER_COMMENT <- "Tutorial comment" #Add info about experiment, eg. explain (Multiple)Comparisons, to include in reports

RUN_PART_TARGETED <- RUN_CODE

#
#####################



##########Targeted_approach##########
if(RUN_PART_TARGETED == RUN_CODE){
  
  ## options
  REIMS <- '.RAW folder from Waters REIMS.'
  EXACTIVE <- '.raw files from Thermo Scientific (Q)-Exactive Hybrid Quadrupole-Orbitrap Mass Spectrometer.'
  
  #Source of file1 variableMetadata:
  VARIABLEMETADATA_FROM_PIPELINE <- 'automatically finds variableMetadata from R pipeline code (based on projectname), present in input folder'
  VARIABLEMETADATA_EXTERN <- 'you choose name of variableMetadata.txt below, present in input folder' #see READ_ME for structure dataframe
  
  VISUALISATION <- "make MZ window plots for each metabolite (width = ppm). Only posible when REIMS pre-processing was performed in R pipeline"
  NO_VISUALISATION <- "make no plots."
  
  CALCULATE_AREA <- "calculate areas. Only posible when REIMS pre-processing was performed in R pipeline"
  NOT_CALCULATE_AREA <- "do not calculate areas."

  
  ## Adjustments
  INSTRUMENT <- REIMS 
  
  #If you choose 'VARIABLEMETADATA_EXTERN', add additional info here:
  VARIABLEMETADATA_EXTERN <- 'TEST_TargetedREIMS_variableMetadata.txt'  #'name.txt' of file. Ignore if file created from pipeline ; name if run modules above: projectname_variableMetadata.txt
  COLLUMN_NR_START_SAMPLES <- 20  #always 20 (auto and manual must be same format); unless extra col merged = 21!
  
  #Source of variableMetadata:
  INPUT_VARIABLES <- VARIABLEMETADATA_EXTERN
  
  #input file3 targetedMetadata:
  INPUT_STANDARDS <- 'targetedMetadata.txt'    #'name.txt'
  PPM <- 100 #e.g. 3 for Q-exactive (lipids), 5 for exacutive and Q-exactive (polar), 50 for REIMS
  
  VISUALISATION_MZ_WINDOW_PPM <- VISUALISATION   #only possible when REIMS pre-processing was performed in R pipeline, since dataset prior to PP is needed.
  AREA <- CALCULATE_AREA  #only possible when REIMS pre-processing was performed in R pipeline, since dataset prior to PP is needed.

}
#
#####################


