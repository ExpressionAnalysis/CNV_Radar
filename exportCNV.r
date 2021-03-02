#!/usr/bin/env Rscript
rm(list=ls())
init <- Sys.time(); timer <- proc.time();

####################################################################################################
# This script is designed to take the standard output .tsv file from CNV Radar and 
#   export the data as a genomic VCF (gVCF)
#
# gVCF conventions are written with the assumption that only one sample per file 
#   is being represented
# Block representation is used with the END of the block reported in the INFO column
#
# No variant sites are represented by a "./." genotype and a "." for the reported alt allele
####################################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Required Packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(getopt) #Consider suppressWarnings(if (!require("getopt")) {install.packages("getopt", dependencies = TRUE)})

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
argString <- commandArgs(trailingOnly = T) # Read in command line arguments

usage <- paste("Usage: Rscript exportCNV.r
               -- Required Parameters --
               [-i | --infile]           <Path to the .tsv output from CNV_Radar.r>
               [-s | --sampleID]         <The name of the sample to put in the VCF header> 
               -- Optional Parameters -- 
               [-f | --format]           <Output format> (Options = 'vcf'; default = 'vcf')
               [-o | --outfile]          <Prefix for the output files> (default = The name of the tsv file up to but not including the .tsv extension)
               [-d | --outpath]          <Path to the directory to save the outputs> (default = current working directory)
               [-p | --purity]           <The tumor purity of the sample, between 0 and 1> (default = 1)
               -- Optional Flags --   
               [-C | --compress]         <gzip the output> (default = F)
               -- Help Flag --  
               [-h | --help]             <Displays this help message>
               Example:
               Rscript exportCNV.r -i <tumor sample>.tsv -s <sample name> -f vcf -o test.out
               \n",sep="")

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(	
  'compress',         'C', 0, "logical",
  'outpath',          'd', 2, "character",
  'infile',           'i', 1, "character",
  'format',           'f', 2, "character",
  'outfile',          'o', 2, "character",
  'purity',           'p', 2, "numeric",
  'sampleID',         's', 1, "character",
  'help',             'h', 0, "logical"
), byrow=TRUE, ncol=4);

#argString <- c("-i", "Z:/data/haven/q805788/HGSC_test_samples/fastqs/TCRBOA6-T-WEX.single_sample_ann.CNVRadar.tsv", "-s", "TCRB06", "-f", "vcf")

params=getopt( spec, argString)

if ( !is.null(params$help) | is.null(params$infile) | is.null(params$sampleID) ) {
  add_to_log(lvl="error", func="getopt", message = "\nEither you asked for help or you are missing a required parameters: infile, sampleID \n\n")
  add_to_log(lvl="error", func="getopt", message = usage)
  q(save="no",status=1,runLast=FALSE)
}

if(is.null(params$outfile)){
  params$outfile = gsub("tsv.*$", "", basename(params$infile), ignore.case = T, perl=T)
} else {
  # If the user supplied output name doesn't end in a . then add it
  if ( substr(x = params$outfile, start = nchar(params$outfile), stop = nchar(params$outfile)) != "." ){ params$outfile <- paste0(params$outfile, ".")}  
}
if(is.null(params$outpath)){params$outpath = getwd()}
if(is.null(params$compress)){params$compress = F}
if(is.null(params$radar_version)){params$radar_version = "1.2.2"} 
if(is.null(params$purity)){params$purity = 1} 
if(is.null(params$format)){params$format = "vcf"}

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
add_to_log(lvl = "info", func="getopt", message=paste0("Params: ", paste(names(params), params, sep=" = ")))

# Ensure that the sampleID is properly formatted
proper_name <- gsub('[^A-Za-z0-9_-]', '', gsub(" ", "_", trimws(params$sampleID)))
if (params$sampleID != proper_name){
  add_to_log(lvl = "info", func="main", message=paste0("Sample ID reformated from ", params$sampleID, " to ", proper_name ))
  params$sampleID <- proper_name
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the required functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
determine_copynumber <- function(x, purity=params$purity, ploidy = 2){
  
  # converts the log2 values to absolute scale
  # The thresholds map to integer copy numbers in order, starting from zero: 
  #     log2 ratios up to the first threshold value are assigned a copy number 0, 
  #     log2 ratios between the first and second threshold values get copy number 1, and so on.
  
  # purity = 1
  # ploidy = 2
  # log2( (1 - purity) + purity * (0:6 + .5) / ploidy )
  
  # Log2 value up to    |   Copy Number
  #       -1.1          |       0       # Empirically set based on the literature
  #       -0.4          |       1
  #        0.3          |       2
  #        0.7          |       3
  #        1.1          |       4
  #        1.4          |       5
  #        1.7          |       6
  #        1.9          |       7
  #        2.1          |       8
  
  # Having more than 6 copies is reported by previous clients as not more functionally relevant than 6, so capping at six, then >6.
  
  tmp <- data.frame("cn" = 0:6, 
                    "threshold" = c(-1.1, log2( (1 - purity) + purity * (1:6 + .5) / ploidy) ), 
                    stringsAsFactors = F)
  
  if (any(x <= tmp$threshold)){
    out <- tmp[min( which( x <= tmp$threshold ) ), "cn"] 
  } else {
    out <- ">6"
  }
  
  return( out )
}

vcf_header <- function(sname = params$sampleID, version = params$radar_version, date = format(Sys.Date(), "%Y%m%d")){
  
  out <- paste0('##fileformat=VCFv4.2
##fileDate=', date, '
##source=CNV_Radar v', version, '
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=FOLD_CHANGE,Number=1,Type=Float,Description="Fold change">
##INFO=<ID=FOLD_CHANGE_LOG,Number=1,Type=Float,Description="Log fold change">
##INFO=<ID=PROBES,Number=1,Type=Integer,Description="Number of probes in CNV">
##ALT=<ID=LOH,Description="Copy Neutral loss of heterozygousity">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t',sname)
  
  return(out)
}

segment2vcf <- function(x){
  # Extract information for the info column
  
  names(x) <- trimws(tolower(names(x)))
  
  # Check to ensure that we have all of the expected columns in the supplied output
  required_cols <- c("log2fc", "iscnv", "isloh_only", "stop", "observeddepth", "qscore", "zscore", "hetvar", "length") 
  if (!all(required_cols %in% names(x) ) ){
    missing_cols <- required_cols[!(required_cols %in% names(x))]
    print(paste0("Supplied CNV Radar output is missing required columns for creating a VCF: ", paste(missing_cols, collapse = ", ")))
    q(save = "no", stats = 1, runLast = F)
  }
  
  log2FC <- as.numeric(x['log2fc'])
  FOLD_CHANGE <- 2^log2FC
  FOLD_CHANGE_LOG <- log2FC
  isCNV <- as.logical(x['iscnv'])
  isLOH <- as.logical(x['isloh_only'])
  END <- x['stop']
  DP <- x['observeddepth']
  QSCORE <- x['qscore']
  ZSCORE <- x['zscore']
  HETVAR <- x['hetvar']
  
  SVLEN <- x['length']
  if (log2FC < 0){SVLEN <- SVLEN * -1}
  
  # Assign the type of CNV
  if ( !isCNV & !isLOH ){
    SVTYPE <- "."
  } else if (isCNV & log2FC < 0){
    SVTYPE <- "DEL"
  } else if (isCNV & log2FC > 0){
    SVTYPE <- "DUP"
  } else if (isLOH) {
    SVTYPE <- "LOH"
  } 
  
  # Determine the absolute copy number 
  #   Assuming 
  #     1. a diploid genome
  #     2. 100% tumor purity
  CN <- determine_copynumber(log2FC)
  
  # Structure the info column
  INFO <- paste("IMPRECISE", 
                paste0("SVTYPE=", SVTYPE),
                paste0("END=", END),
                paste0("SVLEN=", SVLEN),
                paste0("FOLD_CHANGE=", FOLD_CHANGE),
                paste0("FOLD_CHANGE_LOG=", FOLD_CHANGE_LOG),
                paste0("QSCORE=", QSCORE), 
                paste0("ZSCORE=", ZSCORE),
                paste0("HETVAR=", HETVAR),
                sep=";")
  
  # Structure the format field
  FORMAT <- paste("GT", "GQ", sep=":")
  if (SVTYPE == "DUP") {FORMAT <- paste(FORMAT, "CN", "CNQ", sep=":")}
  
  # Structure the sample info field
  if (CN == 0){GT <- "1/1"} else {GT <- "0/1"}
  
  if (SVTYPE == "."){
    GENOTYPE <- paste("./.", QSCORE, sep=":")
  } else if (SVTYPE == "DEL"){
    GENOTYPE <- paste(GT, QSCORE, sep=":")  
  } else if (SVTYPE == "DUP"){
    GENOTYPE <- paste(GT, 0, CN, QSCORE, sep=":")
  } else if (SVTYPE == "LOH"){
    GENOTYPE <- paste(GT, 0, CN, QSCORE, sep=":")
  } 
  
  out <- data.frame("CHROM" = x['chr'],
                    "POS" = x['start'],
                    "ID" = ".",
                    "REF" = "N",
                    "ALT" = SVTYPE,
                    "QUAL" = ".",
                    "FILTER" = ".",
                    "INFO" = INFO,
                    "FORMAT" = FORMAT,
                    "samplename" = GENOTYPE,
                    stringsAsFactors = F)

  return(out)
}

writeVCF <- function(header, body, filename){
  write.table(x=header, file=filename, quote=F, row.names = F, col.names = F, sep = '\t')
  write.table(x=body, file = filename, quote = F, col.names = F, row.names = F, append = T, sep = '\t')  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check that all of the necessary input files exist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!file.exists(params$infile)){
  add_to_log(lvl = "info", func="error", message=paste0("The input file does not exist, please check the supplied path/file extension: ", params$infile))
  q(save = "no", 1, runLast = F)
} 

df <- read.table(params$infile, header=T, sep='\t', stringsAsFactors = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Formatting
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (tolower(params$format) == "vcf"){
  # Alternative methods can be found here https://stackoverflow.com/questions/4227223/convert-a-list-to-a-data-frame
  vcf_body <-  do.call(rbind.data.frame, apply(df, 1, function(x) segment2vcf(x)) )
  colnames(vcf_body)[10] <- params$sampleID
  
  #---------------------------------------------------------------------
  # Basic QC
  #---------------------------------------------------------------------
  if (nrow(df) != nrow(vcf_body)){
    add_to_log(lvl="error", func="segment2vcf", message = paste0("Rows in the VCF (", nrow(df), ") don't match rows in the tsv (", nrow(vcf_body), ")") )
    q(save="no",status=2,runLast=FALSE)
  }
  
  if ( !all(df$Chr == vcf_body$CHROM) | !all(df$Start == vcf_body$POS) ){
    add_to_log(lvl="error", func="segment2vcf", message = paste0("The Chr/Start in the .tsv don't match the CHROM/POS in the .vcf") )
    q(save="no",status=3,runLast=FALSE)
  }
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write the CNV calls out to disk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (tolower(params$format) == "vcf"){
  outfile <- file.path(params$outpath, paste0(params$outfile, "vcf"))
  writeVCF(header=vcf_header(), body=vcf_body, filename=outfile )
  add_to_log(lvl = "info", func="main", message=paste0("VCF file output to ", outfile))
}

if (params$compress & exists("outfile")){
  cmd <- paste("gzip", "-f", outfile) #Use force so that it overwrite any existing files
  system(command = cmd, intern = F, wait = T)
  add_to_log(lvl = "info", func="main", message=paste0("Compressing ", outfile))
} else {
  add_to_log(lvl = "warning", func="main", message="No output file exists to compress")
}

add_to_log(lvl = "info", func="main", message=paste0("Process began at ", init, " and finished at ", Sys.time(), "\n"))
add_to_log(lvl = "info", func="main", message=paste0("Elapsed time: ", (proc.time() - timer)[['elapsed']]))
add_to_log(lvl = "info", func="main", message=("Finished\n"))
