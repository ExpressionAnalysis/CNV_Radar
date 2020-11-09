#!/usr/bin/Rscript
init <- Sys.time(); timer <- proc.time();

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Required Packages
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
                [-i | --indir] < Path to the directory where the roiSummary.txt files are > (Character)
                [-o | --outfile] < Output file name to which '.dendrogram_all.jpeg' will be appended > (Character)
              Optional Parameters:
                [-m | --minD] <Minimum depth required to be considered> (Numeric, Default=10)
                [-p | --outpath] <Path to the directory to save the outputs> (default = current working directory)
                [-r | --roi] < Suffix pattern for the ROI summary file > (Character, Default = '_roiSummary.txt')
              Optional Flags:
                [-h | --help] < Display this help message >
                \n",sep="")

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(
  'outfile',   'o', 1, "character",
  'indir',     'i', 1, "character",
  'minD',      'm', 2, "numeric",
  'outpath',   'p', 2, "character",
  'roi',       'r', 2, "character",
  'help',      'h', 0, "logical"
), byrow=TRUE, ncol=4);

params <- getopt(spec, argString)

# If missing required fields then display usage and quit
if(length(params) < 1 | is.null(params$indir) | is.null(params$outfile)  | !is.null(params$help)) {
  cat("\nYou asked for help or are missing a required parameter: indir, outfile \n\n")
  cat(usage)
  q(save="no",status=1,runLast=FALSE)
}

#Set default values
if(is.null(params$roi)) {params$roi <- "\\_roiSummary\\.txt$"}
if(is.null(params$minD)) {params$minD <- 10}
if(is.null(params$outpath)) {params$outpath <- getwd()}

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the roi summaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(params$indir)
fns <- list.files(pattern = params$roi)

add_to_log(lvl = "debug", func="main", message=paste0(length(fns), " ROI summaries identified in ", params$indir) )

samps <- gsub(pattern = params$roi, replacement = "", x = fns)
d <- NULL
for(j in seq_along(fns)) {
  
  if (j%%5 == 0){
    add_to_log(lvl = "debug", func="main", message=paste0("Reading ROI summary ", j, " of ", length(fns)))
  }
  
  cmd <- paste("cut -f2", fns[j])
  add_to_log(lvl = "debug", func="main", message=paste0("Reading ROI summary ", j, ":", cmd))
  
  tmp <- read.table(file = pipe(cmd), header = T, sep = '\t', as.is = T)
  names(tmp) <- samps[j]
  if(j==1) {
    ids <- read.table(file = pipe(paste("cut -f1", fns[j])), header = T, sep = '\t', as.is = T)
    ids$Chr <- gsub(pattern = "\\_.+", replacement = "", ids$ROI)
    ids$Chr <- gsub(pattern = "^chr", replacement = "", ids$Chr)
    d <- tmp
  } else {
    d <- cbind(d, tmp)
  }
  rm(cmd)
}

add_to_log(lvl = "debug", func="main", message=paste("Completed input of all", length(fns),"ROI summaries"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Subset the data to only those position above the minimum depth 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kp <- apply(X = d, MARGIN = 1, FUN = median) > params$minD & ids$Chr %in% 1:22
d <- d[kp, ]
ids <- ids[kp, ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Median center and log transform the depths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- apply(X = d, MARGIN = 2, FUN = mean)
meds <- median(tmp)
tmp <- matrix(data = tmp, nrow = nrow(d), ncol = ncol(d), byrow = T)
d <- (d / tmp) * meds
d <- log10(d + 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make the dendrogram
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
makeDend <- function(x) {
  dist1 <- 1 - cor(x)
  dist1 <- as.dist(dist1)
  fit <- hclust(dist1, method="ward.D2")
  plot(fit, sub="", main="", xlab="", cex = 0.5)
  return(0)
}

jpeg(filename = file.path(params$outpath, paste0(params$outfile, ".dendrogram_all.jpeg")), width = 8, height = 6, units = "in", quality = 100, res = 300)
makeDend(x = d)
dev.off()

add_to_log(lvl = "info", func="main", message=paste0("Process began at ", init, " and finished at ", Sys.time()))
add_to_log(lvl = "info", func="main", message=paste0("Elapsed time: ", (proc.time() - timer)[['elapsed']]))
add_to_log(lvl = "info", func="main", message=("Finished"))
