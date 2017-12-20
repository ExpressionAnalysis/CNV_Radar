#!/usr/bin/env Rscript

suppressPackageStartupMessages(require("getopt"))
argString <- commandArgs(trailingOnly = TRUE)

usage <- paste("Options:
              Required Parameters:
                [-o | --out] < Path where the output RData file should be saved > (Character)
                [-d | --directory] < Path to the location of the ROI summaries for the normal samples > (Character)
              Optional Parameters:
                [-r | --roi] < Suffix pattern for the ROI summary file > (Character, Default = '_roiSummary.txt')
              Optional Flags:
                [-h | --help] < Display this help message >
                \n",sep="")

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(
  'directory',        'd', 1, "character",
  'out',              'o', 1, "character",
  'roi',       'r', 2, "character",
  'help',             'h', 0, "logical"
), byrow=TRUE, ncol=4);

params <- getopt(spec, argString)

# If missing required fields then display usage and quit
if(length(params) < 1 | is.null(params$out) | is.null(params$directory)  | !is.null(params$help)) {
  cat("\nYou asked for help or are missing a required parameter: out, directory \n\n")
  cat(usage)
  q(save="no",status=1,runLast=FALSE)
}

#Set default values
if(is.null(params$roi)) {params$roi <- "_roiSummary.txt"}

# Read in the depth information file from _roisummary.txt
getROIsum <- function(rois, suffix) {
  d <- NULL # this will be a matrix of ROI X samples
  samps <- gsub(".+/", "", rois)
  samps <- gsub(paste0(suffix,"$"), "", samps)
  # Add another column for each sample
  for(i in 1:length(rois)) {
    print(i)
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

# We can create a control cohort using
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

save(file = params$out, list = c("control", "ids", "RGP"))
print("finished")