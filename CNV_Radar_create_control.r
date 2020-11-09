#!/usr/bin/env Rscript
init <- Sys.time()
timer <- proc.time()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Required packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(require("getopt"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup logging
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_to_log <- function(lvl, func, message){
  # <Date> <function> <level> <information>
  timestamp <- paste0("[",Sys.time(),"]")
  entry <- paste(timestamp, func, toupper(lvl), message, sep = " - ")
  cat(paste0(entry, "\n"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup getopt for command line paramaterization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
argString <- commandArgs(trailingOnly = TRUE)

usage <- paste("Options:
              Required Parameters:
                [-d | --directory] < Path to the location of the ROI summaries for the normal samples > (Character)
              Optional Parameters:
                [-o | --out_file] < Name of the output file, '.RData' will be appended to name provided here > (Character, Default = cnvradar_normal_cohort)
                [-p | --out_dir] < Path where the directory where the RData file should be saved > (Character, Default = getwd())
                [-r | --roi] < Suffix pattern for the ROI summary file > (Character, Default = '_roiSummary.txt')
              Optional Flags:
                [-h | --help] < Display this help message >
                \n",sep="")

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(
  'directory',        'd', 1, "character",
  'out_file',         'o', 2, "character",
  'out_dir',          'p', 2, "character",
  'roi',              'r', 2, "character",
  'help',             'h', 0, "logical"
), byrow=TRUE, ncol=4);

params <- getopt(spec, argString)

# If missing required fields then display usage and quit
if(length(params) < 1 | is.null(params$directory)  | !is.null(params$help)) {
  add_to_log(lvl="warn", func="getopt", message = "\nYou asked for help or are missing a required parameter: directory \n\n")
  add_to_log(lvl="warn", func="getopt", message = usage)
  q(save="no",status=1,runLast=FALSE)
}

#Set default values
if(is.null(params$roi)) {params$roi <- "_roiSummary.txt"}
if(is.null(params$out_file)) {params$out_file <- "cnvradar_normal_cohort"}
if(is.null(params$out_dir)) {params$out_dir <- getwd()}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Required functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getROIsum <- function(rois, suffix) {
  # Read in the depth information file from _roisummary.txt
  d <- NULL # this will be a matrix of ROI X samples
  samps <- gsub(".+/", "", rois)
  samps <- gsub(paste0(suffix,"$"), "", samps)
  samps <- gsub("[_+!@#$?^-]+$", "", samps, perl=T)
  # Add another column for each sample
  for(i in 1:length(rois)) {
    add_to_log(lvl="info", func="getROIsum", message=paste0("Processing control sample:", i ," of ", length(rois)))
    if(is.null(d)) { # if it's the first sample, also include row headers = ROI name
      d <- read.table(file = pipe(paste("cut -f1,2", rois[i])), header = T, sep = "\t", as.is = T)
      names(d)[2] <- samps[i]
    } else {
      tmp <- read.table(file = pipe(paste("cut -f2", rois[i])), header = T, sep = "\t", as.is = T)
      names(tmp) <- samps[i]
      d <- cbind(d, tmp)
    }
  }
  return(d)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Logging information
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
add_to_log(lvl = "info", func="main", message=paste0("User: ", Sys.info()[['effective_user']]) )
add_to_log(lvl = "info", func="main", message=paste0("Running from: ", Sys.info()[['nodename']]) )
add_to_log(lvl = "info", func="main", message=paste0("Platform: ", sessionInfo()$platform) )
add_to_log(lvl = "info", func="main", message=paste0("R version: ", sessionInfo()$R.version$version.string) )
add_to_log(lvl = "info", func="main", message=paste0("R packages loaded: ",  paste(names(sessionInfo()$otherPkgs), collapse=", ") ) )
add_to_log(lvl = "info", func="main", message=paste0("Rscript: ", gsub("--file=", "", grep(pattern = "^--file", commandArgs(trailingOnly = F), value = T))))
add_to_log(lvl = "info", func="main", message=paste0("Arguments: ", paste(commandArgs(trailingOnly = T), collapse=" ")) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a control cohort using supplied samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d <- NULL

if(!is.null(params$directory)) {
  setwd(params$directory)
  rois <- list.files(pattern=params$roi)
  d <- rbind(d, getROIsum(rois, params$roi))
}

# Parse the ROI names into chromosome, beginning and end
colnum <- which(colnames(d) %in% "ROI")
tmp <- gsub(pattern = "^hg[0-9]+\\_", replacement = "", x = d[[colnum]])
tmp <- strsplit(x = tmp, split = "_")
d$Chr <- unlist(lapply(tmp, function(x) { x[1] }))
d$Beg <- as.numeric(unlist(lapply(tmp, function(x) { x[2] })))
d$End <- as.numeric(unlist(lapply(tmp, function(x) { x[3] })))

ids <- d[ ,c("Chr", "Beg", "End")]
if(length(grep("chr",ids$Chr[1]))>0) {
	ids$Chr <- gsub("chr","",ids$Chr)
}

# relative genomic position (RGP) is used for plotting
ids$RGP <- 0
RGP <- list()
for(i in c(1:22, "X", "Y")) {
  RGP[[i]] <- max(ids$RGP)
	ids$RGP[ ids$Chr==i] <- RGP[[i]] + as.numeric(ids$Beg[ ids$Chr==i])
}

# separate the depth data (i.e. "control")
control <- d[ , !names(d) %in% c("Chr", "Beg", "End", "ROI")]
control <- as.matrix(control)

save(file = file.path(params$out_dir, paste0(params$out_file, ".RData")), list = c("control", "ids", "RGP"))
add_to_log(lvl = "info", func="main", message=paste0("Process began at ", init, " and finished at ", Sys.time(), "\n"))
add_to_log(lvl = "info", func="main", message=paste0("Elapsed time: ", (proc.time() - timer)[['elapsed']]))
add_to_log(lvl = "info", func="main", message=("Finished\n"))
