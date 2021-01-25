#!/usr/bin/env Rscript
init <- Sys.time(); timer <- proc.time();

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load Required Packages
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(getopt) #Consider suppressWarnings(if (!require("getopt")) {install.packages("getopt", dependencies = TRUE)})
require(data.table) 
require(yaml) 

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

usage <- paste("Usage: Rscript CNV_Radar.r
               -- Required Parameters --
               [-c | --control]          <Normal control created using CNV_Radar_create_control.r> (Required)
               [-r | --roi]              <Path to the ROI summary> (Required)
               [-v | --vcf]              <Path to the annotated vcf file>(Required)
               -- Optional Parameters -- 
               [-f | --config]           <Path to the yaml formatted config file of model parameters> (default = CNV_Radar_config.yml)
               [-o | --outfile]          <Prefix for the output files> (default = The name of the vcf file up to but not including the .vcf extension)
               [-p | --outpath]          <Path to the directory to save the outputs> (default = current working directory)
               [-d | --depth]            <Comma separated list of format field labels in the VCF file to sum to derive the total read depth> (default = 'DP')
               [-n | --printChrs]        <Comma separated list of Which chromosomes to plot>(default=all)
               [-x | --cex]              <Number indicating the amount by which plotting text and symbols should be scaled relative to the default of 1>(default=0.35)
               -- Optional Flags --   
               [-A | --useAllVars]       <Use all variants instead of only variants with a 'COMMON' annotation from the VCF>(default=FALSE)
               [-G | --gatk]             <Was the variant calling done with GATK> (default = FALSE)
               [-L | --sensitivity]      <Run CNV Radar in standard sensitivity mode>(default=FALSE)
               [-C | --omitCNVcalls]     <Omit lines on the plot to denote called CNVs>(default=FALSE)
               [-D | --plotPredDepth]    <Plot the predicted smoothed log2(FC) depths across the ROI>(default=FALSE)
               [-V | --plotPredVAF]      <Plot the predicted smoothed Heterozygous (VAF-0.5)^2 across the ROI>(default=FALSE)
               -- Output/Reporting Flags --   
               [-W | --writeFilteredVCF] <Write the table of variants after performing filtering>(default=FALSE)					
               -- Help Flag --  
               [-h | --help]             <Displays this help message>
               Example:
               Rscript CNV_Radar.r -r roiSummary.txt -c /data/cnvradar_normal_cohort.RData -o Sample1.cnvradar.out -G
               \n",sep="")

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(	
  'control',          'c', 1, "character",
  'omitCNVcalls',     'C', 0, "logical",
  'depth',            'd', 2, "character",
  'plotPredDepth',    'D', 0, "logical",
  'useAllVars',       'F', 0, "logical",
  'config',           'f', 2, "character",
  'gatk',             'G', 0, "logical", 
  'sensitivity',      'L', 0, "logical",
  'outfile',          'o', 2, "character",
  'outpath',          'p', 2, "character",
  'printChrs',        'n', 2, "character",
  'roi',              'r', 1, "character",
  'plotPredVAF',      'V', 0, "logical",
  'vcf',              'v', 1, "character",
  'writeFilteredVCF', 'W', 0, "logical",
  'cex',              'x', 2, "numeric",
  'help',             'h', 0, "logical"
), byrow=TRUE, ncol=4);

if (length(argString) == 0){
  argString <- c( "-c", "/data/cnvradar_normal_cohort.RData", "-r", "TCRBOA1-T-WEX_roiSummary.txt", "-v", "TCRBOA1-T-WEX.single_sample_ann.vcf.gz", "-o", "TCRBOA1-T-WEX", "-G" ) 
  # "--file", "/scripts/CNV_Radar.r",
}

params=getopt( spec, argString)

if ( !is.null(params$help) | is.null(params$control) | is.null(params$roi) | is.null(params$vcf) ) {
  add_to_log(lvl="error", func="getopt", message = "\nEither you asked for help or you are missing a required parameters: control, roi, vcf\n\n")
  add_to_log(lvl="error", func="getopt", message = usage)
  q(save="no",status=1,runLast=FALSE)
}

# Set default path to the scripts, helper functions and config file
params$wd <- dirname(gsub("--file=", "", grep("--file", commandArgs(trailingOnly = F), value = T)))
params$wd <- gsub("^\\.", getwd(), params$wd) # if the path is relative, get the absolute path 
rm("argString")

# Set defaults for optional processing variables
if(is.null(params$config)){params$config = file.path(params$wd, "CNV_Radar_config.yml") }
if(is.null(params$outfile)){
  params$outfile = gsub("vcf.*$", "", basename(params$vcf), ignore.case = T, perl=T)
} else {
  # If the user supplied output name doesn't end in a . then add it
  if ( substr(x = params$outfile, start = nchar(params$outfile), stop = nchar(params$outfile)) != "." ){ params$outfile <- paste0(params$outfile, ".")}  
}
if(is.null(params$outpath)){params$outpath = getwd()}
if(is.null(params$useAllVars)){params$useAllVars = F}
if(is.null(params$gatk)){params$gatk = F}
if(is.null(params$sensitivity)){params$sensitivity = F}

# Set defaults for optional VCF processing variables
if(is.null(params$depth)){params$depth = "DP"}

# Set defaults for optional plotting variables
if(is.null(params$cex)){params$cex=0.35}
if(is.null(params$printChrs)){
  params$printChrs="all" #default write out the whole genome image
} else {
  # These are the chromosomes that will be plotted as jpeg files    
  if(params$printChrs!="") params$printChrs <- as.numeric(strsplit(x = params$printChrs, split = ",")[[1]]) 
}

# Set defaults for optional NULL plotting parameters
if(is.null(params$omitCNVcalls)){params$omitCNVcalls=F}
if(is.null(params$plotPredDepth)){params$plotPredDepth=F}
if(is.null(params$plotPredVAF)){params$plotPredVAF=F}

# Set defaults for optional reporting variables
if(is.null(params$writeFilteredVCF)){params$writeFilteredVCF=F}

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
# Read in the required functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(file.path(params$wd, "CNV_Radar_functions.r"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check that all of the necessary input files exist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
missing_files <- c()
if (!file.exists(params$config)){missing_files <- c(missing_files, "config")}
if (!file.exists(params$control)){missing_files <- c(missing_files, "control")}
if (!file.exists(params$roi)){missing_files <- c(missing_files, "roi")}
if (!file.exists(params$vcf)){missing_files <- c(missing_files, "vcf")}
if (length(missing_files)!=0){
  add_to_log(lvl = "info", func="error", message=paste0("The following input files do not exist, please check the supplied path/file extension: ", paste(missing_files, collapse = ", ")))
  q(save = "no", 1, runLast = F)
} 
rm(missing_files)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in the modeling config file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model_params <- read_yaml(params$config)

# Choose the overlap threshold to use for the analysis based 
#   on the mode of sensitivity
if(params$sensitivity){
  params$overlap <- model_params$overlap$STANDARD_SENSITIVITY
} else {
  params$overlap <- model_params$overlap$HIGH_SENSITIVITY
}

#--------------------------------------------------------
# Load in the normal control created using 
# CNV_ROI_finder_create_control.r
# This file contains three objects:
# 1. RGP - a list containing a relative base-pair scale,
# 			where the index of the first base on chr2 follows
#				the index for the last base on chr1.  
# 2. ids - a four column dataframe of the regeion of interest
#				1. Chrom ("chr" is stripped from this column)
#				2. Start
#				3. Stop
#				4. Relative Genome Position of the ROI START
# 3. control - an N x M matrix of summarized depths for 
#				N ROI regions and M samples in the normal cohort. 
# 			Should have the same number of rows as item 2 
#--------------------------------------------------------
load(file = params$control) 

# Strip any 'chr' prefix from the ids
ids$Chr <- gsub("^chr", "", ids$Chr, ignore.case = T, perl = T) 

#--------------------------------------------------------
# Read in the VCF file
#--------------------------------------------------------
f_time <- proc.time()
add_to_log(lvl = "debug", func = "readVCF", message = paste("Reading in the vcf file ", params$vcf ))
f <- readVCF(fileName = params$vcf, filter_common = params$useAllVars, depth_fields = params$depth)
add_to_log(lvl = "debug", func = "readVCF", message = paste("VCF file", params$vcf, "with", nrow(f), "variants took", round((proc.time() - f_time)[['elapsed']] / 60,3), "minutes to read and parse"))
rm(f_time)

# Add the relative genomic position for plotting purposes
f$RGP <- determineRGP(f)

if (params$gatk){
  # AF calculations from GATK are biased above max depth (default 250x)
  # Filter out high depth positions
  add_to_log(lvl = "info", func="main", message="Allele frequencies from GATK are biased above 250x depth, removing those positions.")
  if ( !all(is.na(f$DP)) ){
    f <- f[ f$DP < (max(f$DP) * 0.95), ]   
  } else {
    add_to_log(lvl = "warning", func="main", message="Some reported variants were missing depth information, retaining all variants")
  }
}

if(params$writeFilteredVCF) write.table(x = f, file = file.path(params$outpath, paste0(params$outfile, "filtered_vaf.txt")), quote = F, sep = "\t", row.names = F)

#--------------------------------------------------------
# Read in the depths across the region of interest (ROI)
#--------------------------------------------------------
f_time <- proc.time()
add_to_log(lvl = "debug", func = "readROI", message = paste("Reading in the depth file ", params$roi ))
d <- read.table(file = params$roi, header = T, sep = "\t", as.is = T)
add_to_log(lvl = "debug", func = "readROI", message = paste("ROI depth file", params$roi, "with", nrow(d), "regions of interest took", round((proc.time() - f_time)[['elapsed']] / 60,3), "minutes to read and parse"))
rm(f_time)

if (all(colnames(d) != c("ROI", "Mean_Depth"))){
  add_to_log(lvl = "ERROR", func = "readROI", message = paste("The ROI provided is not a 2 column (ROI, Mean_Depth) file. Please check file:", params$roi))
  q(save = no, status = 1, runLast = F)
}

# If the current sample is in the control panel then remove it
#if( grepl(basename(params$sample), colnames(control)) ) {
#  control <- control[ ,-grep(basename(params$sample), colnames(control))]
#}

# Strip any 'chr' prefix from the ROI indicator
d$ROI <- gsub("^chr", "", d$ROI, ignore.case = T, perl = T)

if(nrow(d)==nrow(control)) {
  # Sometimes the ROI positions are in scientific notation in ids
  # Use format(<field>, scientific=F) to ensure identity
  if (all(with(ids, paste(Chr, trimws(format(Beg, scientific=F)),  trimws(format(End, scientific=F)), sep='_')) == gsub("^chr", "", d$ROI, ignore.case = T))){
    # ROI SUMMARY SHOULD ALWAYS HAVE IDENTICAL ROW IDS!!! (else, a merge should be used)
    control <- cbind(d$Mean_Depth, control)   
  } else {
    # FOR THE EXOME, the ROI is too big and the docker runs out of memory and kills the script
    # So you should probably sort the dataframes identically prior to running the script
    add_to_log(lvl = "error", func="ReadROI", message=paste0("The control panel and supplied sample have divergent ROI order, please sort and try again"))
    q(save = "no", status = 1, runLast = F)
    
    # Otherwise this is untested code that I think would work given proper resources
    # Get the names of the columns in the control panel
    control_cols <- colnames(control)
    
    # Setup the merge keys
    control$key <- with(ids, paste(Chr, trimws(format(Beg, scientific=F)),  trimws(format(End, scientific=F)), sep='_'))
    d$ROI <- gsub("^chr", "", d$ROI, ignore.case = T)
    
    # Merge the data
    control <- merge(control, d, by.x = "key", by.y = "ROI", all.x=T)
    
    # Check to see if the merge was successful for every ROI
    if (any(is.na(control$Mean_Depth))){
      add_to_log(lvl = "error", func="ReadROI", message=paste0("The control panel and supplied sample have divergent ROI definitions"))
      q(save = "no", status = 1, runLast = F)
    }
    
    # Reorder/drop columns to prepare for modeling
    control <- control[, c("Mean_Depth", control_cols)]
  }
  
} else {
  add_to_log(lvl = "error", func="ReadROI", message="Control data and roi summary data do not match in number of rows")
  quit(save="no", status=1, runLast=FALSE)
}
rm("d") # general cleanup

control <- control[ order(ids$RGP), ] # Sort according to RGP
ids <- ids[ order(ids$RGP), ] # Sort according to RGP
# save.image("test.RData")

#--------------------------------------------------------
# The variable control is going to be repurposed
# so setup a way to reset it back to the start
#--------------------------------------------------------
depths <- control 

#---------------------------------------------------------------------
# We need to be able to select non-CNV ROIs, 
# at first we don't know where they are so select everything
# suspected to be copy neutral, aka autosomes
#---------------------------------------------------------------------
add_to_log(lvl = "debug", func="main", message="Filtering to autosomes only")
kp <- ids$Chr %in% 1:22 

#---------------------------------------------------------------------
# Iteratively fit the model and adjust parameters accordingly
#---------------------------------------------------------------------
iterative_fit <- list()

# Iteratively increase the sensitivity for identifying CNVs
lowthreshes <- seq(from = model_params$threshLow*model_params$iter, to = model_params$threshLow, length = model_params$iter)
highthreshes <- seq(from = model_params$threshHigh*model_params$iter, to = model_params$threshHigh, length = model_params$iter)


for(i in 1:model_params$iter) {
  iterative_fit$threshLow <- lowthreshes[i] # Set to the ith threshold
  iterative_fit$threshHigh <- highthreshes[i] # Set to the ith threshold
  
  #---------------------------------------------------------------------
  # Median center and log transform the depths
  #---------------------------------------------------------------------
  add_to_log(lvl = "debug", func="main", message=paste("Median normalizing iteration", i))
  control <- depths + 0.5 # to avoid log(0) issues
  medDepths <- apply(control[kp, ], 2, median) 
  control <- control/matrix(data = medDepths, nrow = nrow(control), ncol = ncol(control), byrow = T)
  control <- log(control) 
  
  #---------------------------------------------------------------------
  # Use a linear model to estimate capture bias, using only non-CNV ROIs
  #---------------------------------------------------------------------
  # The first column is the observed depths for the tumor sample (Dependent variable)
  # The other columns are the observed depths from the normal cohort (Independent variable)
  add_to_log(lvl = "debug", func="main", message=paste("Beginning modeling for iteration", i))
  iterative_fit$fit <- lm(formula = control[kp, 1] ~ control[kp, -1]) 
  iterative_fit$fit$fitted.values <- cbind(1, control[, -1]) %*% iterative_fit$fit$coefficients 
  
  #---------------------------------------------------------------------
  # Calculate the log2 fold change
  #   bound the log fold change to avoid highly influential 
  #   (and possibly artifactual) outliers
  #---------------------------------------------------------------------
  add_to_log(lvl = "debug", func="main", message=paste("Calculating fold change for iteration", i))
  ids$fc <- control[,1] -  iterative_fit$fit$fitted.values 
  ids$fc <- pmax(pmin(ids$fc, 3), -3) 
  
  #---------------------------------------------------------------------
  # Initialize empty dataframes for cutpoints and plotting
  #---------------------------------------------------------------------
  cpts <- pall <- NULL
  
  #---------------------------------------------------------------------
  # Analyze each CNV events over each chromosome
  #---------------------------------------------------------------------
  for(chr in unique(ids$Chr)[ grep("^[0-9]", unique(ids$Chr))]) { # Iterate across autosomal chromosomes
    
    add_to_log(lvl = "debug", func="main", message=paste("Analyzing Chr", chr, "for iteration", i))
    if( (sum(ids$Chr==chr) > 10) & (sum(f$CHROM==chr) > 10) ) {
      
      #---------------------------------------------------------------------
      # Going to use a linear spline to predict copy neutrality depth
      #---------------------------------------------------------------------
      # Use a spline to smooth the allele frequency and fold change data
      # Scale the number of knots used to the size of the chromosome
      iterative_fit$cknots <- pmax(20, round(sum(f$CHROM==chr)/model_params$smooth)) 
      
      # Initialize a data frame that will be used as the input to the main function "getCpts"
      preds <- data.frame(Beg=ids$Beg[ ids$Chr==chr], 
                          End=ids$End[ ids$Chr==chr], 
                          RGP=ids$RGP[ ids$Chr==chr], 
                          fc=ids$fc[ ids$Chr==chr], 
                          obs=control[ids$Chr==chr,1], 
                          exp=iterative_fit$fit$fitted.values[ids$Chr==chr], 
                          stringsAsFactors = F)
      
      #---------------------------------------------------------------------
      # When looking at heterozygous variants, what is the deviation
      # away from an allele frequency of 0.5
      # Bound the smoothed values by their range of possible values
      #---------------------------------------------------------------------
      preds$afv <- predict(smooth.spline(x = f$POS[ f$CHROM==chr], 
                                         y = 6*abs(f$AF[ f$CHROM==chr] - 0.5)^3, 
                                         nknots = iterative_fit$cknots, 
                                         w = f$DP[ f$CHROM==chr]), 
                           x = preds$Beg)$y
      preds$afv <- pmax(pmin(preds$afv, 0.75), 0) 
      
      #---------------------------------------------------------------------
      # Determine the smoothed fold change (smfc)
      # Make reasonable bounds to eliminate influential outliers
      #   - deriv is the numerical derivative used to resolve breakpoint focality
      #   - score is the "score" used to identify CNV events
      #---------------------------------------------------------------------
      if( iterative_fit$cknots >= round(nrow(preds)*0.5) ) iterative_fit$cknots <- round(nrow(preds)*0.5) 
      preds$smfc <- predict(smooth.spline(x = preds$Beg, 
                                          y = preds$fc, 
                                          nknots = iterative_fit$cknots, 
                                          w = exp(preds$exp) ), 
                            x = preds$Beg)$y
      preds$smfc <- pmax(pmin(preds$smfc, 3), -3) 
      preds$deriv <- c(diff(preds$smfc), 0) 
      preds$score <- preds$afv * preds$smfc*20 
      
      #---------------------------------------------------------------------
      # Store the predictions for plotting purposes
      #---------------------------------------------------------------------
      if(params$printChrs[1]=="all") pall <- rbind(pall, preds)
      
      #---------------------------------------------------------------------
      # Identify CNV events using the determined score
      #---------------------------------------------------------------------
      chrcpts <- getCpts(preds = preds, chr = chr)
      
      
      #---------------------------------------------------------------------
      # Update kp, so that for the next iteration we are only using ROIs
      #   that are absent of any CNVs based on the current calls
      #---------------------------------------------------------------------
      ###-------------------------------------------###
      ### getCpts returns slightly different tables ###
      ###   depending on high/standard sensitivity  ###
      ###-------------------------------------------###
      if (!params$sensitivity) {
        if(sum(chrcpts$IsCNV) > 0) {
          for(j in which(chrcpts$IsCNV) ) {
            kp[ids$Chr==chrcpts$Chr[j] & chrcpts$Start[j] <= ids$End & ids$Beg <= chrcpts$Stop[j]] <- F
          }
        }
      }else{
        chrcpts <- chrcpts[ , c("Chr", "Start", "Stop", "RGP", "log2FC", "Qscore", "ObservedDepth", "ExpectedDepth", "Zscore", "HetVar")]
        chrcpts <- chrcpts[ chrcpts$ExpectedDepth > model_params$minExpDepth & (chrcpts$Qscore > model_params$threshHigh | chrcpts$Qscore < model_params$threshLow), ] 
        if(nrow(chrcpts) > 0) {
          for(j in 1:nrow(chrcpts)) {
            kp[ids$Chr==chrcpts$Chr[j] & chrcpts$Start[j] <= ids$End & ids$Beg <= chrcpts$Stop[j]] <- F
          }
        }
      }
      
      #---------------------------------------------------------------------
      # Store the checkpoints of the identified CNVs
      #---------------------------------------------------------------------
      cpts <- rbind(cpts, chrcpts)
      
      
      #---------------------------------------------------------------------
      # Make the plot of the log2 fold change by chromosome
      #---------------------------------------------------------------------      
      if(chr %in% params$printChrs & i==model_params$iter) {
        
        add_to_log(lvl = "debug", func = "chrPlot", message = paste("Plotting Chr", chr," for iteration", i))
        
        jpeg(filename = file.path(params$outpath, paste0(params$outfile, "Chr", chr, "_lfc_vaf.jpeg")), width = 8, height = 6, units = "in", quality = 100, res = 300)
        plotCNV(depth = preds, vaf = f, cutpoints = chrcpts)
        dev.off()
        
      }
    }
  }
  
  #---------------------------------------------------------------------
  # For each iteration, plot out a whole genome plot
  #   This is useful for seeing how the threshold for CNVs change with
  #   each iteration
  #---------------------------------------------------------------------
  if(params$printChrs[1]=="all") {
    add_to_log(lvl = "debug", func = "genomePlot", message = paste("Plotting the genome for iteration", i))
    
    
    jpeg(filename = file.path(params$outpath, paste0(params$outfile, "iter", i, "_Genome_lfc_vaf.jpeg")), width = 12, height = 6, units = "in", quality = 100, res = 300)
    plotCNV(depth = pall, vaf = f, cutpoints = cpts)
    dev.off()
    
  }
}

#---------------------------------------------------------------------
# Add in additional information for the output report
#---------------------------------------------------------------------
cpts$RGP <- NULL
cpts$Length <- cpts$Stop - cpts$Start + 1

#---------------------------------------------------------------------
# If run in high sensitivity mode then make LOH determination
#   set true if it is a CNV and low FC, 
#   or if it is not a CNV, but has high HetVar
#---------------------------------------------------------------------
if (!params$sensitivity) {
  cpts$IsLOH_Only <- F 
  cpts$IsLOH_Only[ cpts$IsCNV & model_params$threshFC_Low < cpts$log2FC & cpts$log2FC < model_params$threshFC_High] <- T
  cpts$IsLOH_Only[ !cpts$IsCNV & cpts$HetVar > model_params$LOHthresh & cpts$Length > model_params$LOHminLength ] <- T
}

#---------------------------------------------------------------------
# Make the output table a bit cleaner by rounding
#---------------------------------------------------------------------
for(i in c("log2FC", "Qscore", "ObservedDepth", "ExpectedDepth", "Zscore", "HetVar")) {
  cpts[[i]] <- round(x = cpts[[i]], digits = 3)
}

#---------------------------------------------------------------------
# Write the CNV calls out to disk
#---------------------------------------------------------------------
add_to_log(lvl = "debug", func="main", message=paste("Writing output to:", file.path(params$outpath, paste0(params$outfile, "CNVRadar.tsv")) ))
write.table(x = cpts, file = file.path(params$outpath, paste0(params$outfile, "CNVRadar.tsv")), quote = F, sep = "\t", row.names = F)

add_to_log(lvl = "info", func="main", message=paste0("Process began at ", init, " and finished at ", Sys.time(), "\n"))
add_to_log(lvl = "info", func="main", message=paste0("Elapsed time: ", (proc.time() - timer)[['elapsed']]))
add_to_log(lvl = "info", func="main", message=("Finished\n"))
