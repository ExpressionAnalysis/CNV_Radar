#!/usr/bin/env Rscript

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


# Read in a variant file.  The format changes based on what caller (e.g. vcf from GATK, vcf from VarPROWL, Tute annotation from VarPROWL)
readVCF <- function(params, rgp=NULL) {
  if(params$caller=="Tute") {
    f <- read.table(file = pipe(paste0("zless ", params$sample, params$vcf)), header = T, sep = "\t", as.is = T, quote = )
    f <- f[nchar(f$Ref)==1 & nchar(f$Alt)==1, ] # keep only SNVs
    f <- f[ f$dbSNP!=".", ]
    f <- f[ f$Zygosity=="het", ]
    f$Pos <- f$Start
  } else {
    if(params$caller=="GATK") {
      # Pipe to zless, so we can read in gzipped or uncompressed files with the same command
      f <- read.table(file = pipe(paste0("zless ", params$sample, params$vcf)), header = F, sep = "\t", as.is = T, skip = 1)
      tmp <- lapply(X = f$V10, FUN = function(x) { y <- strsplit(x = x, split = ":")[[1]][2]; y <- strsplit(y, split = ",")[[1]] } )
      f$RefCnts <- as.numeric(unlist(lapply(X = tmp, FUN = function(x) {x[1]}))) # Total reference counts
      f$AltCnts <- as.numeric(unlist(lapply(X = tmp, FUN = function(x) {x[2]}))) # Total variant counts
      f$DP <- f$AltCnts + f$RefCnts # Total depth
      f$AF <- f$AltCnts / f$DP # Allele frequency
      f <- f[ f$DP < (max(f$DP) * 0.95), ] # Filter out high depth positions. AF calculations from GATK are biased above max depth (default 250x)
    } else if(params$caller=="VarPROWL") {
      f <- read.table(file = pipe(paste0("zless ", params$sample, params$vcf)), header = F, sep = "\t", as.is = T)
      f$AF <- suppressWarnings(unlist(lapply(X = f$V10, FUN = function(x) { y <- strsplit(x = x, split = ":")[[1]][2]; y <- as.numeric(y)} )))
      f$DP <- suppressWarnings(unlist(lapply(X = f$V10, FUN = function(x) { y <- strsplit(x = x, split = ":")[[1]][3]; y <- as.numeric(y)} )))
    }
    names(f)[c(1:2, 4:5)] <- c("Chr", "Pos", "Ref", "Alt")
    f$Chr <- gsub("^chr", "", f$Chr)
    f <- f[grep(pattern = ",", x = f$Alt,invert = T), ]  # Grep out complex variants
    f <- f[ nchar(f$Ref)==1 & nchar(f$Alt)==1, ] # Only keep SNPs
    f <- f[ grep("^0[/|]1:", f$V10),] # Keep only positions called as heterozygous
    if(params$filtCommon) f <- f[ unlist(lapply(X = f$V8, FUN = getCommon )) == 1, ]
  }
  # Calculate relative genomic position (RGP), to be used for plotting many chromosomes together
  f$RGP <- f$Pos
  if(!is.null(rgp)) {
    for(i in names(rgp)) {
      f$RGP[ f$Chr==i] <- f$RGP[ f$Chr==i] + rgp[[i]]
    }
  }
  f <- f[ ,c("Chr", "Pos", "DP", "AF", "RGP")]
  return(f)
}

# This function attempts to more accurately resolve breakpoints, based on raw fold change information
# fc = raw fold change; dv = +1 / -1, depending on increase or decrease in CN at the breakpoint from left to right
getFocality <- function(fc, dv, params) {
  mn <- mean(fc)
  if(mn>0) {
    fc <- fc > params$threshHigh
  } else {
    fc <- fc < params$threshLow
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
getCpts <- function(preds, chr) {
  # Use thresholding to get a list of breakpoint candidates, "diff" selects only the first in a series
  cpts <- which(diff(preds$score > params$threshHigh | preds$score < params$threshLow)!=0)
  # We want to resolve the focality of cpts[i] using a window, b-e, around cpts[i]
  if(length(cpts) > 0) {
    for(i in 1:length(cpts)) {
      b <- max(1, cpts[i]-params$focal_thresh)

      if(i>1) b <- max(b, round(cpts[i] - (cpts[i] - cpts[i-1])*0.5) ) #bug fix from v1
      e <- min( nrow(preds), cpts[i]+params$focal_thresh)

        if(i<length(cpts)) e <- min(e, round( cpts[i] + (cpts[i+1]-cpts[i])*0.5)) #bug fix from v1

      if(e-b > 20) {
        cpts[i] <- getFocality(fc = preds$fc[b:e], dv = sign(preds$deriv[cpts[i]]), params = params)+b-1
      }
    }
  }
  cpts <- unique(cpts) #bug fix from v1
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
        if( (mean(rb[1] < a & a < rb[2]) > params$overlap ) |  (mean(ra[1] < b & b < ra[2]) > params$overlap ) ) {
          cpts$e[i] <- cpts$e[i+1]
          cpts <- cpts[ -(i+1), ]
        } else {
          i <- i + 1
        }
        if(i==nrow(cpts)) break
      }
    }
  }
  cpts <- annotateCpts(cpts, preds)
  cpts$IsCNV <- cpts$ExpectedDepth > params$minExpDepth & (cpts$Qscore > params$threshHigh | cpts$Qscore < params$threshLow)
  nonCpts <- cpts[ !cpts$IsCNV , c("b", "e")]
  cpts <- cpts[ cpts$IsCNV , ] # Keep only CNVs passing the threshold
  if(nrow(cpts)>1) {
    cpts$kp <- T
    for(k in 2:nrow(cpts)) {
      if( abs(cpts$log2FC[k-1] - cpts$log2FC[k]) < params$collapse  & ((cpts$b[k] - cpts$e[k-1]) == 1)  ) {
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
  cpts$IsCNV <- cpts$ExpectedDepth > params$minExpDepth & (cpts$Qscore > params$threshHigh | cpts$Qscore < params$threshLow)
  return(cpts)
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
