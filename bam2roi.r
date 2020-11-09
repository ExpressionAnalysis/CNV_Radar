#!/usr/bin/Rscript
rm(list=ls())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Required packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressPackageStartupMessages(require("getopt"))
if (!require("data.table")) {install.packages("data.table", dependencies = TRUE); library("data.table")}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setup getopt for command line parameterization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
argString <- commandArgs(trailingOnly = TRUE)

usage <- paste("Options:
               Required Parameters:
                [-b | --bam] < Path to the bam file for the sample > (Character)
                [-d | --bed] < Path to the bed file describing the capture region of interest > (Character)
               Optional Parameters:
                [-f | --func] < Function to calculate the ROI summary depth >(mean or median, Default='median')
                [-j | --jobsch] < Prefix to use before the command when using a job scheduler> (Character, Default='')
                [-o | --out] < Path to the output directory > (Character, Default=getwd())
               Optional Flags:
                [-z | --verbose] < Print additional criteria for debug purposes >(Default=FALSE)
                [-h | --help] < Display this help message >
                \n",sep="")

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(
  'bam',              'b', 1, "character",
  'bed',              'd', 1, "character",
  'func',             'f', 2, "character",
  'jobsch',           'j', 2, "character",
  'out',              'o', 2, "character",
  'verbose',          'z', 0, "logical",
  'help',             'h', 0, "logical"
), byrow=TRUE, ncol=4);

args=getopt(spec, argString)

# If missing required fields then display usage and quit
if(length(args) < 1 | is.null(args$bam) | is.null(args$bed)  | !is.null(args$help)) {
  cat("\nYou asked for help or are missing a required parameter: bam, bed \n\n")
  cat(usage)
  q(save="no",status=1,runLast=FALSE)
}

#Set default values
if(is.null(args$out)){args$out <- getwd()}
if(is.null(args$verbose)){args$verbose <- FALSE}
if(is.null(args$jobsch)){args$jobsch <- ""}
if(is.null(args$func)){args$func <- "median"}
if(args$verbose){ print(args) }

if(args$func %in% c("mean", "median")){
  cat(paste0("Summarize roi by ", args$func))
} else {
  cat(paste0(args$func, " must be either median or mean"))
  cat(usage)
  q(save="no",status=1,runLast=FALSE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Required functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#-----------------------------------------
# Processing the bam file
#-----------------------------------------

start.timestamp <- timestamp()

bam_name <-  basename(args$bam)

cat(paste0("Processing began at ", start.timestamp, "\n"))
cat(paste0("Generating the ROI summary for the bam file ", bam_name, "\n"))
cat(paste0("The roi is defined by ", args$bed, "\n"))

# The awk statement only holds if you use a three column bedfile
if ( length(strsplit(readLines(file(args$bed, "r"), n=1), "\t")[[1]]) != 3){
  # Cut it to only the first three columns and give it a unique name
  bedfile_cleanup <- T
  cat(paste0("Restricting the bedfile ", args$bed, " to only three columns\n"))
  cmd <- sprintf("cut -f1-3 %s > %s", args$bed, gsub("bam", "three_column.bed", bam_name))
  system(command = cmd, wait=T)
  args$bed <- gsub("bam", "three_column.bed", bam_name)
} else {
  bedfile_cleanup <- F
  cat(paste0("The bedfile ", args$bed, " has only three columns\n"))
}

cmd1 <- sprintf("bedtools intersect -a %s -b %s", args$bam, args$bed)
cmd2 <- sprintf("bedtools genomecov -bga -ibam stdin")
cmd3 <- sprintf("bedtools intersect -wao -a stdin -b %s", args$bed )
cmd4 <- sprintf('awk \'($5!=".") {print $5"_"$6"_"$7,$4,$8}\' /dev/stdin')
outname <- file.path(args$out, gsub(".bam", "_roi.txt.gz", bam_name))

cmd <- sprintf("%s | %s | %s |%s | gzip > %s ", cmd1, cmd2, cmd3, cmd4, outname)
system(command = paste(args$jobsch, cmd), wait = T)

if (bedfile_cleanup){
  cmd <- sprintf("rm %s", args$bed)
  system(command = cmd, wait=F)
}

df <- fread(input = file.path(args$out, paste0(gsub(".bam", "_roi.txt", bam_name),".gz")), sep = ' ', header = F, stringsAsFactors = F)
setnames( df, c("id", "depth", "count") )

if (tolower(args$func) == "mean"){
  #aggregate(test, by="id", FUN=function(x) mean( rep(x["depth"], x["count"]) ), na.rm=T)
  roi <- data.frame("depth"=cbind(by(df, df$id, function(x) mean(rep(x$depth, x$count)))))
} else if (tolower(args$func) == "median"){
  #aggregate(test, by="id", FUN=function(x) median( rep(x["depth"], x["count"]) ), na.rm=T)
  roi <- data.frame("depth"=cbind(by(df, df$id, function(x) median(rep(x$depth, x$count)))))
}

# The CNV_Radar scripts require this header name even though it may be median
colnames(roi) <- c("Mean_Depth")
roi$ROI <- rownames(roi)

write.table(roi[,c("ROI","Mean_Depth")], file.path(args$out, gsub(".bam", "_roiSummary.txt", bam_name)), sep = '\t', col.names = T, row.names = F, quote = F)

cat(paste0("The processing of ", bam_name, " began at ", start.timestamp, "\n"))
cat(paste0("The processing of ", bam_name, " ended at ", timestamp(), "\n"))
cat("Finished")
