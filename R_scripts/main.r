# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Main


##########Configuration##########
#Working directory

## options
PATH_USER <- 'xxx/Pipeline_metabolomics' # example path user
  
CODE_AUTORUN <- 'run code in ternminal automatically'
CODE_DEVELOPMENT <- 'run code manually in Rstudio for development locally'

## Adjustments
PATH <- PATH_USER

CODE_RUN_MODE <- CODE_AUTORUN

#
#####################


##########Sources##########
#Source functions
path_R_scripts <- file.path(PATH, 'R_scripts')
setwd(path_R_scripts)
source('data_loading.R')
source('data_statistics.R')
source('data_converting.R')
source('data_writing.R')
source('data_plot.R')

if (CODE_RUN_MODE == CODE_AUTORUN){
  #Recognize projectname from Rscript commando
  name_project = commandArgs(trailingOnly=TRUE) #name is given in "Rscript main.r projectname" as argument
  # test if there is one argument: if not, return an error
  if (length(args) == 0) {
    stop("Your projectname must be supplied", call.=FALSE)
  } 
  if (length(args) > 1) {
    stop("Only one projectname must be supplied", call.=FALSE)
  }
}
if (CODE_RUN_MODE == CODE_DEVELOPMENT){
  #give projectname in configuration.r (for code in development @Rstudio)
  name_project <- EXPERIMENT
}

#Source configuration and set input folder
path_data_in <- file.path(PATH, 'Data/Input', name_project) #directory must exist!
setwd(path_data_in)
source('Configuration.R')

#make output directory for project:
dir.create(file.path(PATH, 'Data/Output', name_project))
path_data_out <- file.path(PATH, 'Data/Output', name_project)

#
#####################

 
##########R Pipeline - Main##########
print(Sys.time())
start_time_total <- Sys.time()
print("R pipeline - start!")

#Project name:
stopifnot(name_project == EXPERIMENT)
print(name_project)

#Run selected modules
if(RUN_PART_TARGETED == RUN_CODE){
  setwd(path_R_scripts)
  if(INSTRUMENT == REIMS){
    source('Targeted_analysisREIMS.R')
  } else source('Targeted_analysis.R')
}
if(RUN_PART_TARGETED == DONT_RUN_CODE){
  print('no targeted analysis was performed')
}

print(name_project)
print("R pipeline - done!")
print(Sys.time())
end_time_total <- Sys.time()
print(end_time_total - start_time_total)
#
####################
