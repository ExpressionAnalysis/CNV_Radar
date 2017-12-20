#!/usr/bin/env Rscript

require(getopt) #Consider suppressWarnings(if (!require("getopt")) {install.packages("getopt", dependencies = TRUE)})
argString <- commandArgs() # Read in command line arguments

#0=no-arg, 1=required-arg, 2=optional-arg
spec <- matrix(c(	'out',              'O', 1, "character",
					'sample',           'S', 1, "character",
					'threshLow',        'l', 2, "numeric",
					'threshHigh',       'h', 2, "numeric",
					'focal_thresh',     'm', 2, "integer",
					'minVarDepth',      'd', 2, "integer",
					'smooth',           's', 2, "integer",
					'overlap',          'a', 2, "numeric",
					'vcf',              'v', 2, "character",
					'roi',              'R', 1, "character",
					'printChrs',        'P', 2, "character",
					'mmrfReport',       'M', 0, "logical",
					'minExpDepth',      'D', 2, "integer",
					'cex',              'x', 2, "numeric",
					'includeScore',     'u', 0, "logical",
					'iter',             'i', 2, "numeric",
					'filtCommon',       'f', 0, "logical",
					'writeFilteredVCF', 'F', 0, "logical",
					'chr1Stats',        'c', 0, "logical",
					'threshFC_Low',     'L', 2, "numeric",
					'threshFC_High',    'H', 2, "numeric",
					'control',          'C', 1, "character",
					'caller',           'V', 2, "character",
					'collapse',         'p', 2, "numeric",
					'LOHthresh',        'r', 2, "numeric",
					'LOHminLength',     'g', 2, "integer",
					'plotScore',        'e', 0, "logical",
					'plotCNVs',         'n', 0, "logical",
					'highSensitivity',  'I', 0, "logical",
                    'help',             'z', 0, "logical"
), byrow=TRUE, ncol=4);
params=getopt(spec)

if(!is.null(params$help) | is.null(params$out) | is.null(params$sample) | is.null(params$roi) | is.null(params$control)) {
  cat("\nEither you asked for help or you are missing a required parameters\n\n")
  cat(paste("Usage: Rscript CNV_Radar.r
  				-- Required Parameters --
					[-O | --out] <Name of the output file> (Required)
					[-S | --sample] <Name of the sample to process> (Required)
					[-R | --roi] <Extension for the roiSummary> (Required)
					[-C | --control] <Normal control created using CNV_Radar_create_control.r>(Required)
				-- Optional Parameters --
					[-l | --threshLow] <Sensitivity threshold, the predictive score is less than params$threshLow, for detecting cnvs on the last iteration>(default=0.1)
					[-h | --threshHigh] <Sensitivity threshold, Upper threshold for the Q score of a CNV event>(default=0.04)
					[-m | --focal_thresh] <To resolve the focality of the cutpoint using a window of size +/- params$m around the cutpoint>(default=100)
					[-d | --minVarDepth] <Minimum depth required to call a variant>(default=25)
					[-s | --smooth] <Number of knots, max of 20 or the length of the chromosome/params$smooth per chromosome>(default=4)
					[-a | --overlap] <Remove redundant breakpoints, if there is params$overlap in raw fold change estimates between adjacent regions>(default=0.1 {Low Sensitivity} or 0.2 {High Sensitivity})
					[-v | --vcf] <Extenstion of the file name for the vcf file>(default=VarPROWL:'.all_snv.vcf.gz', GATK:'.unifiedGenotyper_5g.snpEFF.vcf.gz')
					[-P | --printChrs] <Which chromosomes to plot>(default=all)
					[-M | --mmrfReport] <MMRF has some specific CNVs associated with it, this produces a report focused on those>(default=FALSE)
					[-D | --minExpDepth] <Minimum predicted depth required to call a CNV event>(default=10)
					[-x | --cex] <Number indicating the amount by which plotting text and symbols should be scaled relative to the default of 1>(default=0.35)
					[-u | --includeScore] <Add the smoothed log2(FC) and Smoothed (VAF-0.5)^2 to the plot>(default=FALSE)
					[-i | --iter] <Number of iterations for reweighting>(default=3)
					[-f | --filtCommon] <Remove common variants from the VCF>(default=TRUE)
					[-F | --writeFilteredVCF] <Write the table of variants after performing filtering>(default=FALSE)
					[-c | --chr1Stats] <Produce a statistical output of CNV events on chromosome 1q>(default=FALSE)
					[-L | --threshFC_Low] <For a called CNV event, to call an LOH event the log2 fold change must be greater than params$theshFC_Low>(default=0.25)
					[-H | --threshFC_High] <For a called CNV event, to call an LOH event the log2 fold change must be greater than params$theshFC_High>(default=0.2)
					[-V | --caller] <What was the caller used to produce the vcf file either VarPROWL or GATK?>(default='GATK')
					[-p | --collapse] <If the difference of the log2FC between two adjacent points are less than params$collapse then merge into one call>(default=0.1)
					[-r | --LOHthresh] <Lower heterozygous allele frequency threshold for calling an LOH event in a 'neutral' region>(default=0.05)
					[-g | --LOHminLength] <Minimum CNV event length threshold in bases for calling an LOH event in a 'neutral' region>(default=2000)
					[-e | --plotScore] <Plot the Score for called CNVs>(default=TRUE)
					[-n | --plotCNVs] <Add a line to the plot to denote called CNVs>(default=TRUE)
					[-I | --highSensitivity] <Run CNV Radar in high sensitivity mode>(default=TRUE)
                    [-z | --help] <Displays this help message>
   Example:
        XXXXXXX
    \n",sep=""))
 q(save="no",status=1,runLast=FALSE)
}

#set defaults for optional NULL args
if(is.null(params$highSensitivity)){params$highSensitivity=T}
if(is.null(params$threshLow)){params$threshLow=-0.1}
if(is.null(params$threshHigh)){params$threshHigh=0.04}
if(is.null(params$focal_thresh)){params$focal_thresh=100}
if(is.null(params$minVarDepth)){params$minVarDepth=25}
if(is.null(params$smooth)){params$smooth=4}
if (params$highSensitivity) {
	if(is.null(params$overlap)){params$overlap=0.2}
}else{
	if(is.null(params$overlap)){params$overlap=0.1}
}
if(is.null(params$vcf)){params$vcf=NULL}
if(is.null(params$printChrs)){params$printChrs="all"}  #default write out the whole genome image
if(params$printChrs!="") params$printChrs <- strsplit(x = params$printChrs, split = ",")[[1]] # These are the chromosomes that will be plotted
if(is.null(params$mmrfReport)){params$mmrfReport=F}
if(is.null(params$minExpDepth)){params$minExpDepth=10}
if(is.null(params$cex)){params$cex=0.35}
if(is.null(params$includeScore)){params$includeScore=F}
if(is.null(params$iter)){params$iter=3}
if(is.null(params$filtCommon)){params$filtCommon=T}
if(is.null(params$writeFilteredVCF)){params$writeFilteredVCF=F}
if(is.null(params$chr1Stats)){params$chr1Stats=F}
if(is.null(params$threshFC_Low)){params$threshFC_Low=-0.25}
if(is.null(params$threshFC_High)){params$threshFC_High=0.2}
if(is.null(params$caller)){params$caller="GATK"}
if(is.null(params$collapse)){params$collapse=0.1}
if(is.null(params$LOHthresh)){params$LOHthresh=0.05}
if(is.null(params$LOHminLength)){params$LOHminLength=2000}
if(is.null(params$plotScore)){params$plotScore=T}
if(is.null(params$plotCNVs)){params$plotCNVs=T}

params$wd <- dirname(substring(argString[grep("--file=", argString)], 8)) # get the base directory of that file
params$wd <- gsub("^\\.", getwd(), params$wd) # if the path is relative, get the absolute path (necessary if the directory changes before sourcing)

source(file.path(params$wd, "CNV_Radar_functions.r")) # Source the helper functions

load(file = params$control) # Load in the normal control created using CNV_ROI_finder_create_control.r

if(!exists("RGP")) RGP <- NULL # Relative genomic position is used for

vcf <- list(VarPROWL=".all_snv.vcf.gz", GATK=".unifiedGenotyper_5g.snpEFF.vcf.gz") # This is the file extension pattern for the vcf files
if(is.null(params$vcf)) params$vcf <- vcf[[params$caller]]
f <- readVCF(params, rgp=RGP) # Read in the vcf file
if(params$writeFilteredVCF) write.table(x = f, file = paste0(params$out, "_vaf.txt"), quote = F, sep = "\t", row.names = F)

#-------------------------------------------------------------------------
# Use bedtools genomecov to make the summary of the region of interest ROI
#-------------------------------------------------------------------------

d <- read.table(file = paste0(params$sample, params$roi), header = T, sep = "\t", as.is = T)
if(length(grep(gsub(".+/", "", params$sample), colnames(control) )) > 0) control <- control[ ,-grep(gsub(".+/", "", params$sample), colnames(control) )]
if(nrow(d)==nrow(control)) {
  control <- cbind(d$Mean_Depth, control) # ROI SUMMARY SHOULD ALWAYS HAVE IDENTICAL ROW IDS!!! (else, a merge should be used)
} else {
  print("Control data and roi summary data do not match in number of rows")
  quit(save="no", status=1, runLast=FALSE)
}
rm(list = c("d", "argString", "vcf")) # Remove unnecessary objects
control <- control[ order(ids$RGP), ] # Sort according to RGP
ids <- ids[ order(ids$RGP), ] # Sort according to RGP

depths <- control # "control" will be median centered (for non-CNV regions) and log transformed
# Iteratively increase the sensitivity for identifying CNVs
lowthreshes <- seq(from = params$threshLow*params$iter, to = params$threshLow, length = params$iter)
highthreshes <- seq(from = params$threshHigh*params$iter, to = params$threshHigh, length = params$iter)

kp <- ids$Chr %in% 1:22 # Initialize the non-CNV ROIs to all autosomal chromosomes
for(i in 1:params$iter) {
  pall <- NULL # This is used for whole genome plots
  params$threshLow <- lowthreshes[i] # Set to the ith threshold
  params$threshHigh <- highthreshes[i] # Set to the ith threshold
  control <- depths + 0.5 # add 0.5 to avoid log(0)
  medDepths <- apply(control[kp, ], 2, median) # Calculate medians among non-CNV ROIs
  control <- control/matrix(data = medDepths, nrow = nrow(control), ncol = ncol(control), byrow = T) # Median center
  control <- log(control) # log transform
  # Use a linear model to estimate capture bias, using only non-CNV ROIs
  params$fit <- lm(formula = control[kp, 1] ~ control[kp, -1]) # The first column is the sample being analyzed - the other columns are the normal cohort
  params$fit$fitted.values <- cbind(1, control[, -1]) %*% params$fit$coefficients # calculate fitted values
  ids$fc <- control[,1] -  params$fit$fitted.values # This is the log-fold change
  ids$fc <- pmax(pmin(ids$fc, 3), -3) # bound the log fold change to avoid highly influential (and possibly artifactual) outliers
  cpts <- NULL
  for(chr in unique(ids$Chr)[ grep("^[0-9]", unique(ids$Chr))]) { # Iterate across autosomal chromosomes
    if( (sum(ids$Chr==chr) > 10) & (sum(f$Chr==chr) > 10) ) {
      # Initialize a data frame that will be used as the input to the main function "getCpts"
      preds <- data.frame(Beg=ids$Beg[ ids$Chr==chr], End=ids$End[ ids$Chr==chr], RGP=ids$RGP[ ids$Chr==chr], fc=ids$fc[ ids$Chr==chr], obs=control[ids$Chr==chr,1], exp=params$fit$fitted.values[ids$Chr==chr], stringsAsFactors = F)
      # Use a spline to smooth the allele frequency and fold change data
      params$cknots <- pmax(20, round(sum(f$Chr==chr)/params$smooth)) # Scale the number of knots used to the size of the chromosome
      # preds$afv is a metric for deviation from 0.5 in variant allele frequency
      preds$afv <- predict(smooth.spline(x = f$Pos[ f$Chr==chr], y = 6*abs(f$AF[ f$Chr==chr] - 0.5)^3, nknots = params$cknots, w = f$DP[ f$Chr==chr]), x = preds$Beg)$y
      preds$afv <- pmax(pmin(preds$afv, 0.75), 0) # Bound the smoothed values by their range of possible values
      # preds$smfc is the smoothed fold change
      if( params$cknots >= round(nrow(preds)*0.5) ) params$cknots <- round(nrow(preds)*0.5) # Scale the number of knots used to the size of the chromosome
      preds$smfc <- predict(smooth.spline(x = preds$Beg, y = preds$fc, nknots = params$cknots, w = exp(preds$exp) ), x = preds$Beg)$y
      preds$smfc <- pmax(pmin(preds$smfc, 3), -3) # Make reasonable bounds to eliminate influential outliers
      preds$deriv <- c(diff(preds$smfc), 0) # preds$deriv is the numerical derivative used to resolve breakpoint focality
      preds$score <- preds$afv * preds$smfc*20 # preds$score is the "score" used to identify CNV events

	  if(params$printChrs[1]=="all") pall <- rbind(pall, preds)
      # getCpts is the main function that uses the score to assign CNV events
      chrcpts <- getCpts(preds = preds, chr = chr)

	  #################
	  ####function getCpts returns slightly different tables depending on high/standard sensitivity
		if (params$highSensitivity) {
			# Update the non-CNV ROIs, based on the current CNV calls
			if(sum(chrcpts$IsCNV) > 0) {
				for(j in which(chrcpts$IsCNV) ) {
					kp[ids$Chr==chrcpts$Chr[j] & chrcpts$Start[j] <= ids$End & ids$Beg <= chrcpts$Stop[j]] <- F
				}
			}
		}else{
			chrcpts <- chrcpts[ , c("Chr", "Start", "Stop", "RGP", "log2FC", "Qscore", "ObservedDepth", "ExpectedDepth", "Zscore", "HetVar")]
			chrcpts <- chrcpts[ chrcpts$ExpectedDepth > params$minExpDepth & (chrcpts$Qscore > params$threshHigh | chrcpts$Qscore < params$threshLow), ] # Keep only CNVs passing the threshold
			# Update the non-CNV ROIs, based on the current CNV calls
			if(nrow(chrcpts) > 0) {
				for(j in 1:nrow(chrcpts)) {
					kp[ids$Chr==chrcpts$Chr[j] & chrcpts$Start[j] <= ids$End & ids$Beg <= chrcpts$Stop[j]] <- F
				}
			}
		}


      cpts <- rbind(cpts, chrcpts)
      # If the current chromosome is in params$printChrs, and we're on the final iteration, make a jpeg plot of the results
      if(chr %in% params$printChrs & i==params$iter) {
        jpeg(filename = paste0(params$out, "_Chr", chr, "_lfc_vaf.jpeg"), width = 8, height = 6, units = "in", quality = 100, res = 300)
        plot(x = preds$Beg, y = log2(exp(preds$fc)), ylim = c(-2.5, 2.5), xlab = "Position", ylab = "Value", main = paste0("Chr", chr, " VAF and ROI depth"), cex=params$cex, pch=16 ) # log of fold change
        abline(h = 0, col = 'gray', lty = 3)
        points(x = f$Pos[ f$Chr==chr & f$DP > params$minVarDepth ], y = f$AF[ f$Chr==chr & f$DP > params$minVarDepth ]+1, col = 'red', cex=params$cex, pch=16) # Allele frequency + 1
        if(params$includeScore) {
          lines(x = preds$Beg, y = log2(exp(preds$smfc)), lwd = 2, col = 'gray')
          points(x = preds$Beg, y = preds$afv*5 - 2, lwd = 2, col = 'orange', cex=params$cex, pch=16)
          legend(x = 15e7, y = -1.5, legend = c("Log2(FC)", "VAF", "Smoothed log2(FC)", "Smoothed (VAF-0.5)^2", "Score (called CNVs)"), fill = c('black', 'red', 'gray', 'orange', 'green'))
        }
        if (params$highSensitivity) {
			if(sum(chrcpts$IsCNV) > 0 & params$plotScore) segments(x0 = chrcpts$Start[chrcpts$IsCNV], x1 = chrcpts$Stop[chrcpts$IsCNV], y0 = chrcpts$log2FC[chrcpts$IsCNV], y1 = chrcpts$log2FC[chrcpts$IsCNV], col = 'green', lwd = 2)
		} else {
			if(nrow(chrcpts) > 0) segments(x0 = chrcpts$Start, x1 = chrcpts$Stop, y0 = chrcpts$log2FC, y1 = chrcpts$log2FC, col = 'green', lwd = 2)
		}
        dev.off()
      }
      if(params$mmrfReport & i==params$iter & chr=="1" & params$chr1Stats) {
        write.table(x = data.frame(Region=c("Chr1_p", "Chr1_q"), rbind(apply(preds[1:12355, c(7, 8, 10)], 2, mean), apply(preds[12356:nrow(preds), c(7, 8, 10)], 2, mean))), file = paste0(params$out, "_1q_stats.tsv"), quote = F, sep = "\t", row.names = F)
      }
    }
  }

  # This is for whole genome plotting
  if(params$printChrs[1]=="all") {
    jpeg(filename = paste0(params$out, "_iter", i, "_Genome_lfc_vaf.jpeg"), width = 12, height = 6, units = "in", quality = 100, res = 300)
    plot(x = pall$RGP, y = log2(exp(pall$fc)), ylim = c(-2.5, 2.5), xlab = "Position", ylab = "Value", main = paste0("Heterozygous allele frequency and log2(fold change)"), cex=params$cex, pch=16, xaxt='n') # log of fold change
    points(x = f$RGP[ f$DP > params$minVarDepth ], y = f$AF[ f$DP > params$minVarDepth ]+1, col = 'red', cex=params$cex, pch=16) # Allele frequency + 1
    if(params$includeScore) {
      lines(x = pall$RGP, y = log2(exp(pall$smfc)), lwd = 2, col = 'gray')
      points(x = pall$RGP, y = pall$afv*5, lwd = 2, col = 'orange', cex=params$cex, pch=16)
      legend(x = 15e7, y = -1.5, legend = c("Log2(FC)", "VAF", "Smoothed log2(FC)", "Smoothed (VAF-0.5)^2", "Score (called CNVs)"), fill = c('black', 'red', 'gray', 'orange', 'green'))
    }
    abline(v = unlist(RGP), col = 'gray', lty = 2)
    abline(h = 0, col = 'darkgray', lty = 3)

	if (params$highSensitivity) {
		if(nrow(cpts) > 0 & params$plotCNVs) segments(x0 = cpts$RGP[cpts$IsCNV], x1 = cpts$Stop[cpts$IsCNV] + cpts$RGP[cpts$IsCNV] - cpts$Start[cpts$IsCNV],
                                                  y0 = cpts$log2FC[cpts$IsCNV], y1 = cpts$log2FC[cpts$IsCNV], col = 'green', lwd = 3)
    } else {
		if(nrow(cpts) > 0) segments(x0 = cpts$RGP, x1 = cpts$Stop + cpts$RGP - cpts$Start, y0 = cpts$log2FC, y1 = cpts$log2FC, col = 'green', lwd = 2)
	}
	dev.off()
  }
}

cpts$RGP <- NULL
cpts$Length <- cpts$Stop - cpts$Start + 1

if (params$highSensitivity) {
	cpts$IsLOH_Only <- F # set true if it is a CNV and low FC, or if it is not a CNV, but has high HetVar
	cpts$IsLOH_Only[ cpts$IsCNV & params$threshFC_Low < cpts$log2FC & cpts$log2FC < params$threshFC_High] <- T
	cpts$IsLOH_Only[ !cpts$IsCNV & cpts$HetVar > params$LOHthresh & cpts$Length > params$LOHminLength ] <- T
}

#########MMRF has some specific CNVs associated with it, this produces a report focused on those
if(params$mmrfReport) {
	if (params$highSensitivity) {
	  for(isLOHonly in c(TRUE, FALSE)) {
		if(isLOHonly) {
		  cptsTMP <- cpts[ cpts$IsLOH_Only==isLOHonly, ]
		} else {
		  cptsTMP <- cpts[ cpts$IsCNV & cpts$IsLOH_Only==isLOHonly, ]
		}
		fn <- paste0(params$out, "_report_IsLOH_Only_", isLOHonly, ".txt")
		# 17p is up to 24000000; P53 is 7668402 - 7687550
		d <- data.frame(PercentChr13Del=0, PercentChr17Del=0, PercentChr17pDel=0, PercentChr1qAmp=0, PercentP53=0, PropVarLess25=mean(f$DP<25), stringsAsFactors = F)
		d$PropROI_Less25=mean(exp(control[,1])*medDepths[1]<25) # QC metric
		d$VAF_Q20 <- quantile(x = f$AF, probs = 0.2) # QC metric
		d$VAF_Q80 <- quantile(x = f$AF, probs = 0.8) # QC metric
		d$VAF_MAD50 <-  median( abs(f$AF-0.5) ) # QC metric
		if(nrow(cptsTMP) > 0) {
		  beg17 <- min(ids$Beg[ ids$Chr==17])
		  end1 <- max(ids$End[ ids$Chr==1])
		  kp <- cptsTMP$log2FC > 0
		  if(isLOHonly) kp <- T
		  d$PercentChr1qAmp <- sum( (pmax(pmin(cptsTMP$Stop, end1), 144013834) - pmin(pmax(cptsTMP$Start, 144013834), end1)) [cptsTMP$Chr==1 & kp]) / (end1 - 144013834)
		  kp <- cptsTMP$log2FC < 0
		  if(isLOHonly) kp <- T
		  d$PercentChr13Del <- sum(cptsTMP$Length[cptsTMP$Chr==13 & kp]) / (max(ids$End[ ids$Chr==13]) - min(ids$Beg[ ids$Chr==13]) + 1)
		  d$PercentChr17Del <- sum(cptsTMP$Length[cptsTMP$Chr==17 & kp]) / (max(ids$End[ ids$Chr==17]) - min(ids$Beg[ ids$Chr==17]) + 1)
		  d$PercentP53 <- sum( (pmax(pmin(cptsTMP$Stop, 7687550), 7668402) - pmin(pmax(cptsTMP$Start, 7668402), 7687550)) [cptsTMP$Chr==17 & kp]) / (7687550 - 7668402)
		  d$PercentChr17pDel <- sum( (pmax(pmin(cptsTMP$Stop, 22066858), beg17) - pmin(pmax(cptsTMP$Start, beg17), 22066858)) [cptsTMP$Chr==17 & kp]) / (22066858 - beg17)
		}
		write.table(x = d, file = fn, quote = F, sep = "\t", row.names = F)
	  }
	} else {
			fn <- paste0(params$out, "_report.txt")
			# 17p is up to 24000000; P53 is 7668402 - 7687550
			d <- data.frame(PercentChr13Del=0, PercentChr17Del=0, PercentChr17pDel=0, PercentChr1qAmp=0, PercentP53=0, PropVarLess25=mean(f$DP<25), stringsAsFactors = F)
			d$PropROI_Less25=mean(exp(control[,1])*medDepths[1]<25) # QC metric
			d$VAF_Q20 <- quantile(x = f$AF, probs = 0.2) # QC metric
			d$VAF_Q80 <- quantile(x = f$AF, probs = 0.8) # QC metric
			d$VAF_MAD50 <-  median( abs(f$AF-0.5) ) # QC metric
		if(nrow(cpts) > 0) {
			beg17 <- min(ids$Beg[ ids$Chr==17])
			end1 <- max(ids$End[ ids$Chr==1])
			d$PercentChr1qAmp <- sum( (pmax(pmin(cpts$Stop, end1), 144013834) - pmin(pmax(cpts$Start, 144013834), end1)) [cpts$Chr==1 & cpts$log2FC > 0]) / (end1 - 144013834)
			d$PercentChr13Del <- sum(cpts$Length[cpts$Chr==13 & cpts$log2FC < 0]) / (max(ids$End[ ids$Chr==13]) - min(ids$Beg[ ids$Chr==13]) + 1)
			d$PercentChr17Del <- sum(cpts$Length[cpts$Chr==17 & cpts$log2FC < 0]) / (max(ids$End[ ids$Chr==17]) - min(ids$Beg[ ids$Chr==17]) + 1)
			d$PercentP53 <- sum( (pmax(pmin(cpts$Stop, 7687550), 7668402) - pmin(pmax(cpts$Start, 7668402), 7687550)) [cpts$Chr==17 & cpts$log2FC < 0]) / (7687550 - 7668402)
			d$PercentChr17pDel <- sum( (pmax(pmin(cpts$Stop, 22066858), beg17) - pmin(pmax(cpts$Start, beg17), 22066858)) [cpts$Chr==17 & cpts$log2FC < 0]) / (22066858 - beg17)
		}
		write.table(x = d, file = fn, quote = F, sep = "\t", row.names = F)
	}

}

for(i in c("log2FC", "Qscore", "ObservedDepth", "ExpectedDepth", "Zscore", "HetVar")) {
  cpts[[i]] <- round(x = cpts[[i]], digits = 3)
}

write.table(x = cpts, file = paste0(params$out, ".tsv"), quote = F, sep = "\t", row.names = F)
print("finished")
