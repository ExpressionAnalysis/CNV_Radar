#!/usr/bin/env Rscript

suppressPackageStartupMessages(require("getopt"))
argString <- commandArgs(trailingOnly = TRUE)

usage <- paste("Options:
              Required Parameters:
                [-b | --bdir] < Path to the directory where the roiSummary.txt files are > (Character)
                [-d | --ds] < Output file name to which '_dendrogram_all.jpeg' will be appended > (Character)
              Optional Parameters:
                [-m | --minD] <Minimum depth required to be considered> (Numeric, Default=10)
                [-r | --roi] < Suffix pattern for the ROI summary file > (Character, Default = '\\_roiSummary\\.txt$')
              Optional Flags:
                [-h | --help] < Display this help message >
                \n",sep="")

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(
  'ds',        'd', 1, "character",
  'bdir',      'b', 1, "character",
  'minD',      'm', 2, "numeric",
  'roi',       'r', 2, "character",
  'help',      'h', 0, "logical"
), byrow=TRUE, ncol=4);

params <- getopt(spec, argString)

# If missing required fields then display usage and quit
if(length(params) < 1 | is.null(params$bdir) | is.null(params$ds)  | !is.null(params$help)) {
  cat("\nYou asked for help or are missing a required parameter: bdir, ds \n\n")
  cat(usage)
  q(save="no",status=1,runLast=FALSE)
}

#Set default values
if(is.null(params$roi)) {params$roi <- "\\_roiSummary\\.txt$"}
if(is.null(params$minD)) {params$minD <- 10}

oldDir <- getwd()
setwd(params$bdir)
fns <- list.files(pattern = params$roi)
samps <- gsub(pattern = params$roi, replacement = "", x = fns)
d <- NULL
for(j in 1:length(fns)) {
  tmp <- read.table(file = pipe(paste("cut -f3", fns[j])), header = T, sep = '\t', as.is = T)
  names(tmp) <- samps[j]
  if(j==1) {
    ids <- read.table(file = pipe(paste("cut -f1", fns[j])), header = T, sep = '\t', as.is = T)
    ids$Chr <- gsub(pattern = "\\_.+", replacement = "", ids$ROI)
    ids$Chr <- gsub(pattern = "^chr", replacement = "", ids$Chr)
    d <- tmp
  } else {
    d <- cbind(d, tmp)
  }
}

kp <- apply(X = d, MARGIN = 1, FUN = median) > params$minD & ids$Chr %in% 1:22
d <- d[kp, ]
ids <- ids[kp, ]
tmp <- apply(X = d, MARGIN = 2, FUN = mean)
meds <- median(tmp)
tmp <- matrix(data = tmp, nrow = nrow(d), ncol = ncol(d), byrow = T)
d <- (d / tmp) * meds
d <- log10(d + 1)

makeDend <- function(x, fn) {
  dist1 <- 1 - cor(x)
  dist1 <- as.dist(dist1)
  fit <- hclust(dist1, method="ward.D2")
  jpeg(filename = fn, width = 8, height = 6, units = "in", quality = 100, res = 300)
  plot(fit, sub="", main="", xlab="", cex = 0.5)
  dev.off()
  return(0)
}

setwd(oldDir)
makeDend(x = d, fn = paste0(params$ds, "_dendrogram_all.jpeg"))
print("finished")