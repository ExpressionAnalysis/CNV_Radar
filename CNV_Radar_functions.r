#!/usr/bin/env Rscript

'%!in%' <- function(x,y)!('%in%'(x,y))

round_any <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

getCommon <- function(x) {
  y <- strsplit(x = x, split = ";")[[1]]
  y <- grep(pattern = "^COMMON", x = y, value = T)
  if(length(y)==0) {
    y <- 0
  } else {
    y <- gsub(pattern = "^COMMON=", replacement = "", x = y)
    y <- max(as.numeric(strsplit(x = y, split = ",")[[1]]))
  }
  return(y)
}

parse_format <- function(x, field){
  tmp <- strsplit(x, split = ":")[[1]]
  f_order <- grep(x=tmp, pattern = paste0("^", field,"$"), ignore.case = F)
  if (length(f_order) == 0){
    f_order <- NA
  }
  return(f_order)
}

parse_sample <- function(x, field){
  if (is.na(field)){
    tmp <- NA
  } else {
    tmp <- strsplit(x, split = ":")[[1]]
    if (length(tmp) >= field){
      tmp <- tmp[field]
    } else {
      add_to_log(lvl="error", func="parse_sample", message = paste0("The field index (", field, ") is greater than the field length (", length(tmp), ")"))
      stop("Error - parse_sample()")
    }
  }
  return(tmp)
}

extract_from_sample <- function(x, field_name, sname){
  if (sname %!in% names(x)){
    add_to_log(lvl="error", func="extract_from_sample", message = paste0("Sample ID ", sname, " not found in VCF."))
    stop("extract_from_sample() - ERROR")
  }
  
  col_index <- parse_format(x = x['FORMAT'], field = field_name)
  z <- parse_sample(x = x[sname], field = col_index)
  return(z)
}

readVCF <- function(fileName, filter_common = params$useAllVars) {
  
  # Check to see that the file exists
  if (!file.exists(fileName)){
    add_to_log(lvl="error", func="readVCF", message = paste0("The VCF ", fileName, " does not exist. Please check the path and try again."))
    q(save = "no", status = 1, runLast = F)
  }
  
  #Read in the VCF one line at a time.  
  con=file(fileName, open="r")
  lin=readLines(con)
  
  #Figure out which lines start with the "#" sign but not the "##" sign
  cNames <- strsplit(as.character(
    substr(lin[setdiff(grep(pattern="#", lin), grep(pattern="##", lin))], 
           start = 2, 
           stop = nchar(lin[setdiff(grep(pattern="#", lin), grep(pattern="##", lin))]))
  ),"\t")
  headers <- lin[grep(pattern="##", lin)]
  
  #Remove the files you don't need anymore
  close(con)
  rm(con, lin)
  
  #Read in the vcf file without reading in the lines starting with #
  vcf <- read.table(fileName, comment.char ="#", sep = '\t', stringsAsFactors = F)
  colnames(vcf) <- cNames[[1]]
  
  add_to_log(lvl = "debug", func = "readVCF", message = paste(nrow(vcf), "total variants reported in the VCF"))
  
  # Remove the 'chr' prefix if present in the CHROM column
  vcf$CHROM <- gsub("^chr", "", vcf$CHROM, ignore.case = T, perl = T) 
  
  # Keep only the variants on autosomal or sex chromosomes
  autosomes <- vcf$CHROM %in% c(seq(1,22,1),"X", "Y")
  add_to_log(lvl = "debug", func = "readVCF", message = paste(sum(autosomes), "variants reported on the autosomal or sex chromosomes"))
  vcf <- vcf[autosomes, ]
  
  # Filter out complex variants
  complex_vars <- grep(pattern = ",", x = vcf$ALT)
  add_to_log(lvl = "debug", func = "readVCF", message = paste(length(complex_vars), "complex variants to remove from VCF"))
  vcf <- vcf[-complex_vars, ]
  
  # Filter out INDELs and retain only SNVs
  only_snps <- nchar(vcf$REF)==1 & nchar(vcf$ALT)==1
  add_to_log(lvl = "debug", func = "readVCF", message = paste(sum(!only_snps), "indels to remove from VCF"))
  vcf <- vcf[ only_snps, ]
  
  # Keep only reported heterozygous positions
  vcf$GT <- apply(vcf, 1, function(x) extract_from_sample(x, field_name = "GT", sname = colnames(vcf)[10] ))
  het_genotypes <- c("0|1", "0/1", "1|0", "1/0") # include both phased (|) and unphased (/) genotypes
  hom_genotypes <- c("1|1", "1/1")
  add_to_log(lvl = "debug", func = "readVCF", message = paste(sum(vcf$GT %in% het_genotypes), "heterozygous and", sum(vcf$GT %in% hom_genotypes), "homozygous variants reported"))
  vcf <- vcf[vcf$GT %in% het_genotypes,]
  
  # Extract out the depth at each variant position
  vcf$DP <- as.numeric(apply(vcf, 1, function(x) extract_from_sample(x, field_name = "DP", sname = colnames(vcf)[10] )))
  
  # Extract out the allele frequency
  vcf$AF <- as.numeric(apply(vcf, 1, function(x) extract_from_sample(x, field_name = "AF", sname = colnames(vcf)[10] )))
  
  if ( all(is.na(vcf$AF)) ){
    # The allele frequency wasn't there, attempting to calculate it using depths
    # It is assumed that the AD column exists and lists the alleleic depths for the ref (first) and alt (last) alleles
    #     separated by a comma
    vcf$AD <- apply(vcf, 1, function(x) extract_from_sample(x, field_name = "AD", sname = colnames(vcf)[10]))
    vcf$RefCnts <- as.numeric(sapply(strsplit(vcf$AD,","), `[`, 1))
    vcf$AltCnts <- as.numeric(sapply(strsplit(vcf$AD,","), `[`, 2)) 
    vcf$DP <- vcf$RefCnts + vcf$AltCnts
    vcf$AF <- vcf$AltCnts / vcf$DP
  }
  
  # Should we further filter down to variants annotated to be 'COMMON'
  if (!filter_common){
    add_to_log(lvl="debug", func="readVCF", message = "Filtering VCF to only 'common' variants")
    vcf$iscommon <- sapply(X = vcf$INFO, FUN = getCommon ) == 1
    add_to_log(lvl = "debug", func = "readVCF", message = paste(sum(vcf$iscommon), "variants annotated as common"))
    vcf <- vcf[vcf$iscommon, ]
  } 
  
  add_to_log(lvl = "debug", func = "readVCF", message = paste(nrow(vcf), "total heterozygous common single nucleotide variants reported in the VCF"))
  
  chr_breakdown <- table(factor(vcf$CHROM, levels=c(1:22,"X", "Y")))
  add_to_log(lvl = "debug", func = "readVCF", message = paste0("Variants reported on chromosome ", paste(names(chr_breakdown), chr_breakdown, sep=" = ")))
  
  return(vcf[, c(colnames(vcf)[1:10], "GT","AF", "DP")])
}

determineRGP <- function(tmp, rgp=NULL){
  
  if ( any(c("CHROM", "POS") %!in% colnames(tmp))){
    add_to_log(lvl="error", func="determineRGP", message = paste0(paste(c("CHROM", "POS")[c("CHROM", "POS") %!in% colnames(tmp)], collapse = ", "), " not a column in VCF."))
    stop("determineRGP() - ERROR")
  }
  
  # Calculate relative genomic position (RGP), to be used for plotting many chromosomes together
  tmp$RGP <- tmp$POS
  
  if(!is.null(RGP)) {
    for(i in names(RGP)) {
      tmp$RGP[ tmp$CHROM==i] <- tmp$RGP[ tmp$CHROM==i ] + RGP[[i]]
    }
  }
  
  return(tmp$RGP)
}

# This function attempts to more accurately resolve breakpoints, based on raw fold change information
# fc = raw fold change; dv = +1 / -1, depending on increase or decrease in CN at the breakpoint from left to right
getFocality <- function(fc, dv, threshLow, threshHigh) {
  mn <- mean(fc)
  if(mn>0) {
    fc <- fc > threshHigh
  } else {
    fc <- fc < threshLow
  }
  m1 <- cumsum(fc) / 1:length(fc)
  m2 <- rev(cumsum(rev(fc)) / 1:length(fc))
  if(sign(mn)==dv) {
    ans <- m2 - m1
  } else {
    ans <- m1 - m2
  }
  ans[c(1:10, (length(ans): (length(ans) - 9) ))] <- 0
  return(which.max( ans ))
}

# This is the main function that calculates where the breakpoints are, resolves the focality of breakpoints,
# removes redundant breakpoints, and assigns significance of the regions defined by the breakpoints by
# comparison with a normal cohort.
getCpts <- function(preds, chr, 
                    iter_threshLow = iterative_fit$threshLow, 
                    iter_threshHigh = iterative_fit$threshHigh, 
                    overlap=params$overlap,
                    minExpDepth = model_params$minExpDepth,
                    focal_thresh = model_params$focal_thresh,
                    collapse = model_params$collapse,
                    sensitivity = params$sensitivity) {
  # Use thresholding to get a list of breakpoint candidates, "diff" selects only the first in a series
  cpts <- which(diff(preds$score > iter_threshHigh | preds$score < iter_threshLow)!=0)
  
  # We want to resolve the focality of cpts[i] using a window, b-e, around cpts[i]
  if(length(cpts) > 0) {
    for(i in 1:length(cpts)) {
      b <- max(1, cpts[i]-focal_thresh)
      if(i>1) b <- max(b, round((cpts[i] - cpts[i-1])*0.5) )
      e <- min( nrow(preds), cpts[i]+focal_thresh)
      if(i<length(cpts)) e <- min(e, round( (cpts[i+1]-cpts[i])*0.5))
      if(e-b > 20) {
        cpts[i] <- getFocality(fc = preds$fc[b:e], 
                               dv = sign(preds$deriv[cpts[i]]), 
                               threshLow = iter_threshLow,
                               threshHigh = iter_threshHigh)+b-1
      }
    }
  }
  
  if(length(cpts)==0) cpts <- nrow(preds) # if no cpts found
  
  if(cpts[length(cpts)]!=nrow(preds) )  cpts <- c(cpts, nrow(preds) ) # add the last position, if not there
  # remove redundant breakpoints, by seeing if there is a lot of overlap in raw fold change estimates between
  # adjacent regions
  cpts <- data.frame(b=c(1, cpts[1:length(cpts)-1]+1), e=cpts, stringsAsFactors=F) # change to begin, end matrix
  for(k in 1:2) {
    if(nrow(cpts) > 1) {
      i <- 1
      while(TRUE) {
        a <- preds$fc[cpts$b[i]:cpts$e[i]]
        b <- preds$fc[cpts$b[i+1]:cpts$e[i+1]]
        ra <- quantile(x = a, probs = c(0.25, 0.75))
        rb <- quantile(x = b, probs = c(0.25, 0.75))
        if( (mean(rb[1] < a & a < rb[2]) > overlap ) |  (mean(ra[1] < b & b < ra[2]) > overlap ) ) {
          cpts$e[i] <- cpts$e[i+1]
          cpts <- cpts[ -(i+1), ]
        } else {
          i <- i + 1
        }
        if(i==nrow(cpts)) break
      }
    }
  }
  
  if (sensitivity){
    cpts <- cpts[cpts$e - cpts$b > 1, ] # Remove CNVs that have length 1
    # Annotate the regions defined by the breakpoints with log2FC, ObservedDepth, Expected Depth, etc.
    cpts$Qscore <- cpts$Start <- cpts$Stop <- cpts$RGP <- cpts$log2FC <- cpts$ObservedDepth <- cpts$ExpectedDepth <- cpts$Zscore <- cpts$SE <- cpts$HetVar <- NA
    cpts$Chr <- chr
    for(k in 1:nrow(cpts)) { # k <- 1
      cpts$HetVar[k] <- mean(preds$afv[ cpts$b[k]:cpts$e[k]])
      cpts$Start[k] <- preds$Beg[ cpts$b[k]]
      cpts$Stop[k] <- preds$End[ cpts$e[k]]
      cpts$RGP[k] <- preds$RGP[ cpts$b[k]]
      cpts$ExpectedDepth[k] <- mean(exp(preds$exp[cpts$b[k]:cpts$e[k]]) * medDepths[1] )
      cpts$ObservedDepth[k] <- mean(exp(preds$obs[cpts$b[k]:cpts$e[k]]) * medDepths[1] )
      cpts$log2FC[k] <- log2(exp(median(preds$fc[cpts$b[k]:cpts$e[k]]) ))
      all_fc_k <- apply(control[ids$Chr==chr, ][cpts$b[k]:cpts$e[k], ], 2, mean)
      cpts$SE[k] <- sd(all_fc_k[-1])
      cpts$Zscore[k] <- ( all_fc_k[1] - mean(all_fc_k[-1]) ) / cpts$SE[k]
      cpts$Qscore[k] <- mean(preds$score[ cpts$b[k]:cpts$e[k]])
    }
    return(cpts)
  } else {
    cpts <- annotateCpts(cpts, preds)
    
    cpts$IsCNV <- cpts$ExpectedDepth > minExpDepth & (cpts$Qscore > iter_threshHigh | cpts$Qscore < iter_threshLow)
    nonCpts <- cpts[ !cpts$IsCNV , c("b", "e")]
    cpts <- cpts[ cpts$IsCNV , ] # Keep only CNVs passing the threshold
    if(nrow(cpts)>1) {
      cpts$kp <- T
      for(k in 2:nrow(cpts)) {
        if( abs(cpts$log2FC[k-1] - cpts$log2FC[k]) < collapse  & ((cpts$b[k] - cpts$e[k-1]) == 1)  ) {
          cpts$b[k] <- cpts$b[k-1]
          cpts$kp[k-1] <- F
        }
      }
      cpts <- cpts[ cpts$kp, ]
      cpts$kp <- NULL
    }
    cpts <- rbind(cpts[ ,c("b", "e")], nonCpts[ ,c("b", "e")])
    cpts <- cpts[ order(cpts$b), ]
    cpts <- annotateCpts(cpts, preds)
    cpts <- cpts[ , c("Chr", "Start", "Stop", "RGP", "log2FC", "Qscore", "ObservedDepth", "ExpectedDepth", "Zscore", "HetVar")]
    cpts$IsCNV <- cpts$ExpectedDepth > minExpDepth & (cpts$Qscore > iter_threshHigh | cpts$Qscore < iter_threshLow)
    return(cpts)
  }
}

# Annotate the regions defined by the breakpoints with log2FC, ObservedDepth, Expected Depth, etc.
annotateCpts <- function(cpts, preds) {
  cpts$Qscore <- cpts$Start <- cpts$Stop <- cpts$RGP <- cpts$log2FC <- cpts$ObservedDepth <- cpts$ExpectedDepth <- cpts$Zscore <- cpts$SE <- cpts$HetVar <- NA
  cpts$Chr <- chr
  for(k in 1:nrow(cpts)) {
    cpts$HetVar[k] <- mean(preds$afv[ cpts$b[k]:cpts$e[k]])
    cpts$Start[k] <- preds$Beg[ cpts$b[k]]
    cpts$Stop[k] <- preds$End[ cpts$e[k]]
    cpts$RGP[k] <- preds$RGP[ cpts$b[k]]
    cpts$ExpectedDepth[k] <- mean(exp(preds$exp[cpts$b[k]:cpts$e[k]]) * medDepths[1] )
    cpts$ObservedDepth[k] <- mean(exp(preds$obs[cpts$b[k]:cpts$e[k]]) * medDepths[1] )
    cpts$log2FC[k] <- log2(exp(median(preds$fc[cpts$b[k]:cpts$e[k]]) ))
    if(cpts$e[k] - cpts$b[k] > 0) {
      all_fc_k <- apply(control[ids$Chr==chr, ][cpts$b[k]:cpts$e[k], ], 2, mean)
    } else {
      all_fc_k <- control[ids$Chr==chr, ][cpts$b[k]:cpts$e[k], ]
    }
    cpts$SE[k] <- sd(all_fc_k[-1])
    cpts$Zscore[k] <- ( all_fc_k[1] - mean(all_fc_k[-1]) ) / cpts$SE[k]
    cpts$Qscore[k] <- median(preds$score[ cpts$b[k]:cpts$e[k]])
  }
  return(cpts)
}

# Function to add chromosomme annotations to the whole genome plot
find_midpoint <- function(x, end){
  # This function takes a vector of values and finds the midpoint between each index
  
  midpoint <- rep(NA, length(x))
  for (i in seq_along(x)){
    if (i != length(x)){
      midpoint[i] <- median(c(x[i], x[i+1]))
    } else {
      # For the very last row
      midpoint[i] <- median(c(x[i], end))
    }
  }
  
  return(midpoint)
}

even_breaks <- function(max_pnt, brkpnts = 10){
  brks <- seq(0, max_pnt, by = 10^(floor(log10(max_pnt / brkpnts))))
  
  # Find how we can get evenly spaced intervals
  no_remainder <- c()
  for ( i in seq_along(brks)){
    if (length(brks) %% i == 0){no_remainder <- append(no_remainder, i)}
  }
  
  if (length(no_remainder) > 2){
    brks <- brks[seq(1, length(brks), by = length(brks) / max(head(no_remainder, -1)))]  
  } else {
    brks <- brks[seq(1, length(brks), by = length(brks) / tail(no_remainder, 1))]  
  }
  
  return(brks)
}

plotCNV <- function(depth, 
                    vaf, 
                    cutpoints, 
                    relGenomPos = RGP,
                    chr = params$printChrs, 
                    cex_size = params$cex, 
                    depthThresh = model_params$minVarDepth, 
                    omit_cnv_score = params$omitCNVcalls,
                    include_smfc = params$plotPredDepth,
                    include_smvaf = params$plotPredVAF){
  
  if (tolower(chr) == "all"){
    xcol_dp <- xcol_vaf <- xcol_cpts <- "RGP"
    plot_title_prefix <- "Genome wide"
    
    # Subset the variant allele frequency to exclude sex chromosomes
    vaf <- vaf[toupper(vaf$CHROM) %!in% c("X", "Y"), ]
  } else {
    xcol_dp <- "Beg"
    xcol_vaf <- "POS"
    xcol_cpts <- "Start"
    plot_title_prefix <- paste0("Chr", chr)
    
    # Subset variant allele frequency to the plotting focus
    vaf <- vaf[vaf$CHROM == chr, ]
  }
  
  # Initialize the legend placeholders.
  legend_label <- legend_color <- c()
  
  # Plot the log2 fold change of the sequencing depth
  plot(x = depth[, xcol_dp], 
       y = log2(exp(depth$fc)), 
       ylim = c(-2.5, 2.5), 
       xlab = "", 
       ylab = "", 
       main = paste(plot_title_prefix, "heterozygous allele frequency and log2(fold change)"), 
       cex = cex_size, 
       pch = 16, 
       las = 1,
       xaxt='n') 
  
  # Add in the axis label for the black dots
  title(xlab="Position", line=0, cex.lab=0.8, col.lab = "black")
  title(ylab="Depth\nLog(Fold Change)", line=2, cex.lab=0.8, col.lab = "black")
  
  legend_label <- append(x = legend_label, values = "Log2(FC)")
  legend_color <- append(x = legend_color, values = "black")
  
  # Plot the allele frequency but shift it up 1 so not to clutter the visualization
  points(x = vaf[vaf$DP > depthThresh, xcol_vaf], 
         y = vaf[vaf$DP > depthThresh, "AF"]+1, 
         col = 'red', 
         cex = cex_size, 
         pch = 16) 
  
  # Add in the y axis label for the red dots
  title(ylab=paste0(paste(rep(" ", 53), collapse=""), "VAF + 1"), line=2, cex.lab=0.8, col.lab = "red")
  
  legend_label <- append(x = legend_label, values = "VAF")
  legend_color <- append(x = legend_color, values = "red")
  
  # Add in the threshold lines
  if(include_smfc) {
    lines(x = depth[,xcol_dp], 
          y = log2(exp(depth$smfc)), 
          lwd = 2, 
          col = 'gray')
    legend_label <- append(x = legend_label, values = "Smoothed log2(FC)")
    legend_color <- append(x = legend_color, values = "gray")
  }
  
  if(include_smvaf) {
    points(x = depth[,xcol_dp], 
           y = depth$afv*5, 
           lwd = 2, 
           col = 'orange', 
           cex=cex_size, 
           pch=16)
    
    legend_label <- append(x = legend_label, values = "Smoothed (VAF-0.5)^2")
    legend_color <- append(x = legend_color, values = "orange")
  }
  
  # Add in a line to identify the segments where a CNV was identified
  if (!params$sensitivity) {
    # High Sensitivity
    if(sum(cutpoints$IsCNV) > 0 & !omit_cnv_score) {
      segments(x0 = cutpoints[cutpoints$IsCNV, xcol_cpts], 
               x1 = cutpoints[cutpoints$IsCNV, xcol_cpts] + cutpoints$Stop[cutpoints$IsCNV] - cutpoints$Start[cutpoints$IsCNV],
               y0 = cutpoints$log2FC[cutpoints$IsCNV], 
               y1 = cutpoints$log2FC[cutpoints$IsCNV], 
               col = 'green', 
               lwd = 3)
    }
  } else {
    # Standard Sensitivity
    if(nrow(cutpoints) > 0 & !omit_cnv_score) {
      segments(x0 = cpts[,xcol_cpts], 
               x1 = cpts[,xcol_cpts] + cpts$Stop - cpts$Start, 
               y0 = cpts$log2FC, 
               y1 = cpts$log2FC, 
               col = 'green', 
               lwd = 2)
    } 
  }
  
  legend_label <- append(x = legend_label, values = "Score (called CNVs)")
  legend_color <- append(x = legend_color, values = "green")
  
  # Add the line for a copy neutral depth fold change
  abline(h = 0, col = 'darkgray', lty = 3)
  
  if (tolower(chr) == "all"){
    # Add lines to show the division between chromosomes
    abline(v = unlist(relGenomPos), col = 'gray', lty = 2)
  
    # Add Chromosome number to the plot for genome plot
    text(x = find_midpoint(x = unlist(relGenomPos), end = max(c(vaf$RGP, depth$RGP)))[1:22], 
                           y = 2.35, 
                           labels = names(relGenomPos)[1:22], 
                           col = "blue",
                           cex = 0.75)
  } else {
    # Add uniformly spaced lines across the plot
    gridlines <- even_breaks(max(c(vaf[,xcol_vaf], depth[,xcol_dp])))
    
    abline(v = gridlines, col = 'gray', lty = 2)
    
    # Label the number of bases
    text(x = gridlines,
         y = 2.35, 
         labels = formatC(gridlines, format = "e", digits = 1), 
         col = "blue",
         cex = 0.5)
  }
    
  # Add a legend
  legend(x = 15e7, 
         y = -1.5, 
         legend = legend_label, 
         fill = legend_color)
  
}
