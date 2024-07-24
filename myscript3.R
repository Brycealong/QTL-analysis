# #!usr/bin/env Rscript
# Args <- commandArgs(trailingOnly = TRUE)
# 
# if (length(Args) < 3){
#   stop("At least three argument must be supplied (input vcf) (high bulk name) (low bulk name)", call. = FALSE)
# } else if (length(Args) == 3){
#   # default output file
#   message("using the default output file name 'out.png' ...")
#   Args[4] <- "out.png"
# }
# # store the current directory
# initial.dir <- getwd()
# # change to the new directory
# setwd("/home/bryce/Desktop/QTL_analysis/")

# load the necessary libraries
suppressPackageStartupMessages({library(dplyr, quietly = T)
library(tidyr, quietly = T)
library(ggplot2, quietly = T)
library(optparse, quietly = T)
# library(patchwork)
library(vcfR, quietly = T)})


#load the package
#library(QTLseqr)

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "input VCF file", metavar = "character"),
  make_option(c("-h", "--highBulk"), type = "character", default = NULL, 
              help = "high bulk name", metavar = "character"),
  make_option(c("-l", "--lowBulk"), type = "character", default = NULL, 
              help = "low bulk name", metavar = "character"),
  make_option(c("--removeChromList"), type = "character", default = NULL, 
              help = "a list of chromosomes to be removed in the analysis, separated by \",\"", metavar = "character"),
  make_option(c("--refAlleleFreq"), type = "numeric", default = 0.2, 
              help = "filter by reference allele frequency: [refAlleleFreq] <= REF_PRQ <= 1-[refAlleleFreq] [default %default]", metavar = "numeric"),
  make_option(c("--filterAroundMedianDepth"), type = "numeric", default = NULL, 
              help = "filter around median depth", metavar = "numeric"),
  make_option(c("--minTotalDepth"), type = "numeric", default = 100, 
              help = "filter by total sample read depth >= [minTotalDepth] [default %default]", metavar = "numeric"),
  make_option(c("--maxTotalDepth"), type = "numeric", default = 400, 
              help = "filter by total sample read depth <= [maxTotalDepth] [default %default]", metavar = "numeric"),
  make_option(c("--minSampleDepth"), type = "numeric", default = 40, 
              help = "filter by per sample read depth >= [minSampleDepth] [default %default]", metavar = "numeric"),
  make_option(c("--depthDifference"), type = "numeric", default = 100, 
              help = "filter by depth difference <= [depthDifference] between high and low bulks [default %default]", metavar = "numeric"),
  make_option(c("--minGQ"), type = "numeric", default = 99, 
              help = "filter by Genotype Quality: GQ >= [minGQ] [default %default]", metavar = "numeric"),
  make_option(c("--windowSize"), type = "numeric", default = 1e6, 
              help = "window size for QTLseq analysis [default %default]", metavar = "numeric"),
  make_option(c("--popStruc"), type = "character", default = "F2", 
              help = "population structure for QTLseq analysis, either \"F2\" or \"RIL\"", metavar = "character"),
  make_option(c("--bulkSize"), type = "character", default = "25", 
              help = "bulk size(s) for QTLseq analysis. If the bulks are of different sizes(e.g. 385 and 430 for high and low bulk respectively), we set --bulkSize=385,430. If your bulks are the same size you can simply set one value, i.e. --bulkSize=25.", metavar = "character"),
  # make_option(c("--depth"), type = "character", default = NULL, 
  #             help = "min and max depths for QTLseq simulation, --depth=min,max. If not defined, use min and max depth from data", metavar = "character"),
  make_option(c("--replications"), type = "numeric", default = 10000, 
              help = "number of replications for QTLseq simulation [default %default]", metavar = "numeric"),
  make_option(c("--filter"), type = "numeric", default = 0.3, 
              help = "Keeping SNPs with >= [filter] SNP-index in both simulated bulks [default %default]", metavar = "numeric"),
  make_option(c("--intervals"), type = "character", default = "95,99", 
              help = "confidence intervals for QTLseq analysis, [default %default]", metavar = "character"),
  make_option(c("--subset"), type = "character", default = NULL, 
              help = "subset of chromosomes to plot, separated by \",\"", metavar = "character"),
  make_option(c("--plotVar"), type = "character", default = "deltaSNP", 
              help = "variable to plot, either \"deltaSNP\" or \"nSNPs\"", metavar = "character"),
  make_option(c("--plotIntervals"), action = "store_true", default = FALSE, 
              help = "plot confidence intervals"),
  make_option(c("--outputPlot"), type = "character", default = "my_BSA_QTL.png", 
              help = "output plot file name [default %default]", metavar = "character"),
  make_option(c("--exportCSV"), action = "store_true", default = FALSE, 
              help = "export QTL table to CSV"),
  make_option("--intervalCSV", type = "numeric", default = 95,
              help = "threshold for exporting the QTL table. Should match one of the intervals calculated above. [default %default]", metavar = "numeric"),
  make_option(c("--outputCSV"), type = "character", default = "my_BSA_QTL.csv", 
              help = "output CSV file name for QTL table [default %default]", metavar = "character")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$input) || is.null(opt$highBulk) || is.null(opt$lowBulk)) {
  print_help(opt_parser)
  stop("At least three arguments must be supplied (input VCF file, high bulk name, low bulk name).", call. = FALSE)
}

#Set sample and file names
file <- opt$input
highBulk <- unlist(strsplit(opt$highBulk, ","))
lowBulk <- unlist(strsplit(opt$lowBulk, ","))
removeChromList <- if (!is.null(opt$removeChromList)) unlist(strsplit(opt$removeChromList, ",")) else NULL

# Filtering parameters
refAlleleFreq <- opt$refAlleleFreq
filterAroundMedianDepth <- opt$filterAroundMedianDepth
minTotalDepth <- opt$minTotalDepth
maxTotalDepth <- opt$maxTotalDepth
minSampleDepth <- opt$minSampleDepth
depthDifference <- opt$depthDifference
minGQ <- opt$minGQ

# QTLseq analysis parameters
windowSize <- opt$windowSize
popStruc <- opt$popStruc
bulkSize <- as.numeric(unlist(strsplit(opt$bulkSize, ",")))
# depth <- if (!is.null(opt$depth)) seq(as.numeric(unlist(strsplit(opt$depth, ",")[[1]]))) else NULL
replications <- opt$replications
filter <- opt$filter
intervals <- as.numeric(unlist(strsplit(opt$intervals, ",")))

# Plotting parameters
subset <- if (!is.null(opt$subset)) unlist(strsplit(opt$subset, ",")) else NULL
plotVar <- opt$plotVar
plotIntervals <- opt$plotIntervals
outputFig <- opt$outputPlot

# QTL table export parameters
exportCSV <- opt$exportCSV
intervalCSV <- opt$intervalCSV
outputCSV <- opt$outputCSV


mergeBulk <- function(df, bulk) {
  df %>% 
    dplyr::filter(Indiv %in% bulk) %>% 
    dplyr::select(-Indiv) %>%
    tidyr::separate(
      col = "AD",
      into = c("AD_REF", "AD_ALT"),
      sep = ",",
      # extra = "merge",
      convert = TRUE
    ) %>%
    group_by(Key) %>%
    summarise(across(c(AD_REF, AD_ALT, DP), sum),
              across(c(GQ), mean))
}

importFromVCF <- function(file,
                          highBulk,
                          lowBulk,
                          removeChromList = NULL) {
  
  vcf <- vcfR::read.vcfR(file = file)
  # message("Keeping SNPs that pass all filters")
  # vcf <- vcf[vcf@fix[, "FILTER"] == "PASS"] 
  message("Keeping SNPs that are biallelic")
  vcf <- vcf[vcfR::is.biallelic(vcf),]
  
  fix <- dplyr::as_tibble(vcf@fix[, c("CHROM", "POS", "REF", "ALT")]) %>% dplyr::mutate(Key = seq(1:nrow(.)))
  
  # if (!all(
  #     c(
  #         "CHROM", 
  #         "POS", 
  #         paste0(highBulk, ".AD"), 
  #         paste0(lowBulk, ".AD"), 
  #         paste0(highBulk, ".DP"), 
  #         paste0(lowBulk, ".DP")
  #     ) %in% names(SNPset))) {
  #     stop("One of the required fields is missing. Check your VCF file.")
  # }
  
  tidy_gt <- extract_gt_tidy(vcf, format_fields = c("AD", "DP", "GQ"), gt_column_prepend = "", alleles = FALSE)
  
  lowBulkSNPset <- mergeBulk(tidy_gt, lowBulk)
  highBulkSNPset <- mergeBulk(tidy_gt, highBulk)
  SNPset <- lowBulkSNPset %>%
    dplyr::left_join(highBulkSNPset,
                     by = "Key",
                     suffix = c(".LOW", ".HIGH")) %>%
    dplyr::full_join(x = fix, by = "Key") %>%
    dplyr::mutate(
      POS = as.numeric(POS),
      AD_ALT.HIGH = DP.HIGH - AD_REF.HIGH,
      AD_ALT.LOW = DP.LOW - AD_REF.LOW,
      SNPindex.HIGH = AD_ALT.HIGH / DP.HIGH,
      SNPindex.LOW = AD_ALT.LOW / DP.LOW,
      REF_FRQ = (AD_REF.HIGH + AD_REF.LOW) / (DP.HIGH + DP.LOW),
      deltaSNP = SNPindex.HIGH - SNPindex.LOW
    ) %>%
    dplyr::select(-Key)
  #Keep only wanted chromosomes
  if (!is.null(removeChromList)) {
    message("Removing the following chromosomes: ", paste(removeChromList, collapse = ", "))
    SNPset <- SNPset[!SNPset$CHROM %in% removeChromList, ]
  }
  na.omit(as.data.frame(SNPset))
}

filterSNPs <- function(SNPset,
                       refAlleleFreq,
                       filterAroundMedianDepth,
                       minTotalDepth,
                       maxTotalDepth,
                       minSampleDepth,
                       depthDifference,
                       minGQ,
                       verbose = TRUE) {
  
  org_count <- nrow(SNPset)
  count <- nrow(SNPset)
  
  # Filter by total reference allele frequency
  if (!is.null(refAlleleFreq)) {
    if (verbose) {
      message(
        "Filtering by reference allele frequency: ",
        refAlleleFreq,
        " <= REF_FRQ <= ",
        1 - refAlleleFreq
      )
    }
    SNPset <- dplyr::filter(SNPset, SNPset$REF_FRQ < 1 - refAlleleFreq &
                              SNPset$REF_FRQ > refAlleleFreq)
    if (verbose) {
      message("...Filtered ", count - nrow(SNPset), " SNPs")
    }
    count <- nrow(SNPset)
  }
  
  #Total read depth filtering
  
  if (!is.null(filterAroundMedianDepth)) {
    # filter by Read depth for each SNP FilterByMAD MADs around the median
    madDP <-
      mad(
        x = (SNPset$DP.HIGH + SNPset$DP.LOW),
        constant = 1,
        na.rm = TRUE
      )
    medianDP <-
      median(x = (SNPset$DP.HIGH + SNPset$DP.LOW),
             na.rm = TRUE)
    maxDP <- medianDP + filterAroundMedianDepth * madDP
    minDP <- medianDP - filterAroundMedianDepth * madDP
    SNPset <- dplyr::filter(SNPset, (DP.HIGH + DP.LOW) <= maxDP &
                              (DP.HIGH + DP.LOW) >= minDP)
    
    if(verbose) {message("Filtering by total read depth: ",
                         filterAroundMedianDepth,
                         " MADs arround the median: ", minDP, " <= Total DP <= ", maxDP)
      message("...Filtered ", count - nrow(SNPset), " SNPs")}
    count <- nrow(SNPset)
    
  }
  
  if (!is.null(minTotalDepth)) {
    # Filter by minimum total SNP depth
    if (verbose) {
      message("Filtering by total sample read depth: Total DP >= ",
              minTotalDepth)
    }
    SNPset <-
      dplyr::filter(SNPset, (DP.HIGH + DP.LOW) >= minTotalDepth)
    
    if (verbose) {
      message("...Filtered ", count - nrow(SNPset), " SNPs")
    }
    count <- nrow(SNPset)
  }
  
  if (!is.null(maxTotalDepth)) {
    # Filter by maximum total SNP depth
    if (verbose) {
      message("Filtering by total sample read depth: Total DP <= ",
              maxTotalDepth)
    }
    SNPset <-
      dplyr::filter(SNPset, (DP.HIGH + DP.LOW) <= maxTotalDepth)
    if (verbose) {
      message("...Filtered ", count - nrow(SNPset), " SNPs")
    }
    count <- nrow(SNPset)
  }
  
  
  # Filter by min read depth in either sample
  if (!is.null(minSampleDepth)) {
    if (verbose) {
      message("Filtering by per sample read depth: DP >= ",
              minSampleDepth)
    }
    SNPset <-
      dplyr::filter(SNPset, DP.HIGH >= minSampleDepth &
                      DP.LOW >= minSampleDepth)
    if (verbose) {
      message("...Filtered ", count - nrow(SNPset), " SNPs")
    }
    count <- nrow(SNPset)
  }
  
  # Filter by Genotype Quality
  if (!is.null(minGQ)) {
    if (all(c("GQ.LOW", "GQ.HIGH") %in% names(SNPset))) {
      if (verbose) {
        message("Filtering by Genotype Quality: GQ >= ", minGQ)
      }
      SNPset <-
        dplyr::filter(SNPset, GQ.LOW >= minGQ & GQ.HIGH >= minGQ)
      if (verbose) {
        message("...Filtered ", count - nrow(SNPset), " SNPs")
      }
      count <- nrow(SNPset)} 
    else {
      message("GQ columns not found. Skipping...")
    }
  }
  
  # Filter by difference between high and low bulks
  if (!is.null(depthDifference)) {
    if (verbose) {
      message("Filtering by difference between bulks <= ",
              depthDifference)
    }
    SNPset <-
      dplyr::filter(SNPset, abs(DP.HIGH - DP.LOW) <= depthDifference)
    if (verbose) {
      message("...Filtered ", count - nrow(SNPset), " SNPs")
    }
    count <- nrow(SNPset)
  }
  
  if (verbose) {
    message(
      "Original SNP number: ",
      org_count,
      ", Filtered: ",
      org_count - count,
      ", Remaining: ",
      count
    )
  }
  return(as.data.frame(SNPset))
}


countSNPs <- function(POS, windowSize) {
  nout <- length(POS)
  left <- 1
  right <- 1
  out <- numeric(nout)
  
  for (i in 1:nout) {
    while ((right < nout) && (POS[right + 1] <= POS[i] + windowSize / 2)) {
      right <- right + 1
    }
    
    while (POS[left] <= POS[i] - windowSize / 2) {
      left <- left + 1
    }
    
    out[i] <- right - left + 1
  }
  
  return(out)
}

tricubeStat <- function(POS, Stat, windowSize = 2e6, ...)
{
  if (windowSize <= 0)
    stop("A positive smoothing window is required")
  stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0), ...), POS)
}

# tricubeStat_2 <- function(POS, Stat, windowSize = 2e6, windowStepSize = 1e4, ...) {
#   # Check if the windowSize is positive
#   if (windowSize <= 0)
#     stop("A positive smoothing window is required")
#   
#   # Fit the local regression model using locfit
#   fit <- locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0), ...)
#   
#   # Define the prediction points based on step_size
#   prediction_points <- seq(min(POS), max(POS), by = windowStepSize)
#   
#   # Generate predictions at specified points
#   predictions <- stats::predict(fit, prediction_points)
#   
#   # Return the predictions
#   return(predictions)
# }


simulateAlleleFreq <- function(n, pop = "F2") {
  if (pop == "F2") {
    mean(sample(
      x = c(0, 0.5, 1),
      size = n,
      prob = c(1, 2, 1),
      replace = TRUE
    ))
  } else {
    mean(sample(
      x = c(0, 1),
      size = n,
      prob = c(1, 1),
      replace = TRUE
    ))
  }
}

simulateSNPindex <-
  function(depth,
           altFreq1,
           altFreq2,
           replicates = 10000,
           filter = NULL) {
    
    SNPindex_H <- rbinom(replicates, size = depth, altFreq1) / depth
    SNPindex_L <- rbinom(replicates, size = depth, altFreq2) / depth
    deltaSNP <- SNPindex_H - SNPindex_L
    
    if (!is.null(filter)) {
      deltaSNP <- deltaSNP[SNPindex_H >= filter | SNPindex_L >= filter]
    }
    deltaSNP
  }

simulateConfInt <- function(popStruc = "F2",
                            bulkSize,
                            depth = 1:100,
                            replications = 10000,
                            filter = 0.3,
                            intervals = c(0.05, 0.025)) {
  if (popStruc == "F2") {
    message(
      "Assuming bulks selected from F2 population, with ",
      paste(bulkSize, collapse = " and "),
      " individuals per bulk."
    )
  } else {
    message(
      "Assuming bulks selected from RIL population, with ",
      bulkSize,
      " individuals per bulk."
    )
  }
  
  if (length(bulkSize) == 1) {
    message("The 'bulkSize' argument is of length 1, setting number of individuals in both bulks to: ", bulkSize)
    bulkSize[2] <- bulkSize[1]
  }
  
  if (length(bulkSize) > 2) {
    message("The 'bulkSize' argument is larger than 2. Using the first two values as the bulk size.")
  }
  
  if (any(bulkSize < 0)) {
    stop("Negative bulkSize values")
  }
  
  #makes a vector of possible alt allele frequencies once. this is then sampled for each replicate
  tmp_freq <-
    replicate(n = replications * 10, simulateAlleleFreq(n = bulkSize[1], pop = popStruc))
  
  tmp_freq2 <-
    replicate(n = replications * 10, simulateAlleleFreq(n = bulkSize[2], pop = popStruc))
  
  message(
    paste0(
      "Simulating ",
      replications,
      " SNPs with reads at each depth: ",
      min(depth),
      "-",
      max(depth)
    )
  )
  message(paste0(
    "Keeping SNPs with >= ",
    filter,
    " SNP-index in both simulated bulks"
  ))
  
  # tmp allele freqs are sampled to produce 'replicate' numbers of probablities. these
  # are then used as altFreq probs to simulate SNP index values, per bulk.
  CI <- sapply(
    X = depth, 
    FUN = function(x)
    {
      quantile(
        x = simulateSNPindex(
          depth = x,
          altFreq1 = sample(
            x = tmp_freq,
            size = replications,
            replace = TRUE
          ),
          altFreq2 = sample(
            x = tmp_freq2,
            size = replications,
            replace = TRUE
          ),
          replicates = replications,
          filter = filter
        ),
        probs = intervals,
        names = TRUE
      )
    }
  )
  
  CI <- as.data.frame(CI)
  
  if (length(CI) > 1) {
    CI <- data.frame(t(CI))
  }
  
  names(CI) <- paste0("CI_", 100 - (intervals * 200))
  CI <- cbind(depth, CI)
  
  #to long format for easy plotting
  # tidyr::gather(data = CI,
  #     key = interval,
  #     convert = TRUE,
  #     value = SNPindex,-depth) %>%
  #     dplyr::mutate(Confidence = factor(ifelse(
  #         interval > 0.5,
  #         paste0(round((1 - interval) * 200, digits = 1), "%"),
  #         paste0((interval * 200), "%")
  # )))
  CI
}

runQTLseqAnalysis <- function(SNPset,
                              windowSize = 1e6,
                              popStruc = "F2",
                              bulkSize,
                              depth = NULL,
                              replications = 10000,
                              filter = 0.3,
                              intervals = c(95, 99),
                              ...) {
  
  message("Counting SNPs in each window...")
  SNPset <- SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(nSNPs = countSNPs(POS = POS, windowSize = windowSize))
  
  message("Calculating tricube smoothed delta SNP index...")
  SNPset <- SNPset %>%
    dplyr::mutate(tricubeDeltaSNP = tricubeStat(POS = POS, Stat = deltaSNP, windowSize, ...))
  
  #convert intervals to quantiles
  if (all(intervals >= 1)) {
    message(
      "Returning the following two sided confidence intervals: ",
      paste(intervals, collapse = ", ")
    )
    quantiles <- (100 - intervals) / 200
  } else {
    stop(
      "Convidence intervals ('intervals' paramater) should be supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95% confidence interval, 2.5% on each side."
    )
  }
  
  #calculate min depth per snp between bulks
  SNPset <-
    SNPset %>%
    dplyr::mutate(minDP = pmin(DP.LOW, DP.HIGH))
  
  SNPset <-
    SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(tricubeDP = floor(tricubeStat(POS, minDP, windowSize = windowSize, ...)))
  
  if (is.null(depth)) {
    message(
      "Variable 'depth' not defined, using min and max depth from data: ",
      min(SNPset$minDP),
      "-",
      max(SNPset$minDP)
    )
    depth <- min(SNPset$minDP):max(SNPset$minDP)
  }
  
  #simualte confidence intervals
  CI <-
    simulateConfInt(
      popStruc = popStruc,
      bulkSize = bulkSize,
      depth = depth,
      replications = replications,
      filter = filter,
      intervals = quantiles
    )
  
  
  #match name of column for easier joining of repeat columns
  names(CI)[1] <- "tricubeDP"
  
  #use join as a quick way to match min depth to matching conf intervals.
  SNPset <-
    dplyr::left_join(x = SNPset,
                     y = CI #, commented out becuase of above change. need to remove eventually
                     # by = c("tricubeDP" = "depth"))
    )
  as.data.frame(SNPset)
  
}

format_genomic <- function(...) {
  # Format a vector of numeric values according
  # to the International System of Units.
  # http://en.wikipedia.org/wiki/SI_prefix
  #
  # Based on code by Ben Tupper
  # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
  # Args:
  #   ...: Args passed to format()
  #
  # Returns:
  #   A function to format a vector of strings using
  #   SI prefix notation
  #
  
  function(x) {
    limits <- c(1e0, 1e3, 1e6)
    #prefix <- c("","Kb","Mb")
    
    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)
    
    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)
    
    paste(format(round(x/limits[i], 1),
                 trim=TRUE, scientific=FALSE, ...)
          #  ,prefix[i]
    )
  }
}

format_prefix <- 
  function(x) {
    limits <- c(1e0, 1e3, 1e6)
    prefix <- c("","Kb","Mb")
    
    # Vector with array indices according to position in intervals
    i <- findInterval(abs(x), limits)
    
    # Set prefix to " " for very small values < 1e-24
    i <- ifelse(i==0, which(limits == 1e0), i)
    prefix[i]
  }

plotQTLStats <-
  function(SNPset,
           subset = NULL,
           var = "deltaSNP",
           line = TRUE,
           plotIntervals = FALSE,
           ...) {
    
    if (!all(subset %in% unique(SNPset$CHROM))) {
      whichnot <-
        paste(subset[base::which(!subset %in% unique(SNPset$CHROM))], collapse = ', ')
      stop(paste0("The following are not true chromosome names: ", whichnot))
    }
    
    if (!var %in% c("nSNPs", "deltaSNP"))
      stop(
        "Please choose one of the following variables to plot: \"nSNPs\", \"deltaSNP\""
      )
    
    SNPset <-
      if (is.null(subset)) {
        SNPset
      } else {
        SNPset[SNPset$CHROM %in% subset,]
      }
    
    p <- ggplot2::ggplot(data = SNPset) +
      ggplot2::scale_x_continuous(breaks = seq(from = 0, to = max(SNPset$POS), by = 10^(floor(log10(max(SNPset$POS))))), 
                                  labels = format_genomic(), 
                                  name = paste0("Genomic Position (", format_prefix(max(SNPset$POS)), ")")) +
      ggplot2::theme(plot.margin = ggplot2::margin(
        b = 10,
        l = 20,
        r = 20,
        unit = "pt"
      ))
    
    if (var == "nSNPs") {
      p <- p + ggplot2::ylab("Number of SNPs in window")
    }
    
    if (var == "deltaSNP") {
      var <- "tricubeDeltaSNP"
      p <-
        p + ggplot2::ylab(expression(Delta * " (SNP-index)")) +
        ggplot2::ylim(-0.55, 0.55) +
        ggplot2::geom_hline(yintercept = 0,
                            color = "black",
                            alpha = 0.4)
      if (plotIntervals == TRUE) {
        
        ints_df <-
          dplyr::select(SNPset, CHROM, POS, dplyr::matches("CI_")) %>% 
          tidyr::pivot_longer(
            cols = dplyr::matches("CI_"),
            names_to = "Interval", 
            values_to = "value"
          )
        
        p <- p + 
          ggplot2::geom_line(data = ints_df, ggplot2::aes(
            x = POS, 
            y = value, 
            color = Interval
          )) +
          ggplot2::geom_line(data = ints_df, ggplot2::aes(
            x = POS,
            y = -value,
            color = Interval
          ))
      }
    }
    
    if (line) {
      p <-
        p + ggplot2::geom_line(ggplot2::aes(.data[["POS"]], .data[[var]]), ...)
    }
    
    if (!line) {
      p <-
        p + ggplot2::geom_point(ggplot2::aes(.data[["POS"]], .data[[var]]), ...)
    }
    
    p <- p + ggplot2::facet_wrap(~ CHROM, scales = "free_x", ncol = 5)
    p
  }

getQTLTable <-
  function(SNPset,
           interval = 95,
           export = FALSE,
           fileName = "QTL.csv"){
    conf <- paste0("CI_", interval)
    
    if (any(names(SNPset) %in% conf)) {
      stop(
        "Cant find the requested confidence interval. Please check that the requested interval exsits or first use runQTLseqAnalysis to calculate confidence intervals"
      )
    }
    
    #QTL <- getSigRegions(SNPset = SNPset, method = method, interval = interval, alpha = alpha)
    #merged_QTL <- dplyr::bind_rows(QTL, .id = "id")
    SNPset <- SNPset %>%
      dplyr::group_by(CHROM)
    
    qtltable <-
      SNPset %>% 
      dplyr::mutate(passThresh = abs(tricubeDeltaSNP) > abs(!!as.name(conf))) %>%
      dplyr::group_by(CHROM, run = {
        run = rle(passThresh)
        rep(seq_along(run$lengths), run$lengths)
      }) %>%
      dplyr::filter(passThresh == TRUE) %>% 
      dplyr::ungroup() %>%
      dplyr::group_by(CHROM) %>% 
      dplyr::group_by(CHROM, qtl = {
        qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
      }) %>%
      #dont need run variable anymore
      dplyr::select(-run) %>%
      dplyr::summarize(
        start = min(POS),
        end = max(POS),
        length = end - start,
        nSNPs = length(POS),
        avgSNPs_Mb = round(length(POS) / (max(POS) - min(POS)) * 1e6),
        peakDeltaSNP = ifelse(
          mean(tricubeDeltaSNP) >= 0,
          max(tricubeDeltaSNP),
          min(tricubeDeltaSNP)
        ),
        posPeakDeltaSNP = POS[which.max(abs(tricubeDeltaSNP))],
        avgDeltaSNP = mean(tricubeDeltaSNP)
      )
    
    qtltable <- as.data.frame(qtltable)
    
    if (export) {
      write.csv(file = fileName,
                x = qtltable,
                row.names = FALSE)
    }
    return(qtltable)
  }

#Import SNP data from file
df <-
  importFromVCF(
    file,
    highBulk,
    lowBulk,
    removeChromList
  )

# ggplot(data = df) + 
#   geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
#   xlim(0,1000)
# 
# ggplot(data = df) +
#   geom_histogram(aes(x = REF_FRQ))

#Filter SNPs based on some criteria
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = refAlleleFreq,
    filterAroundMedianDepth = filterAroundMedianDepth,
    minTotalDepth = minTotalDepth,
    maxTotalDepth = maxTotalDepth,
    minSampleDepth = minSampleDepth,
    depthDifference = depthDifference,
    minGQ = minGQ,
    verbose = TRUE
  )

# ggplot(data = df_filt) + 
#   geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
#   xlim(0,1000)
# 
# ggplot(data = df_filt) +
#   geom_histogram(aes(x = REF_FRQ))

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = windowSize,
  popStruc = popStruc,
  bulkSize = bulkSize,
  depth = NULL,
  replications = replications,
  filter = filter,
  intervals = intervals
)

#Plot
p <- plotQTLStats(SNPset = df_filt, var = plotVar, plotIntervals = plotIntervals)

ggplot2::ggsave(filename = outputFig, plot = p, width = 8, height = 5)

#export summary CSV
if (exportCSV){
  getQTLTable(SNPset = df_filt, interval = intervalCSV, export = TRUE, fileName = outputCSV)
}
# change back to the original directory
# setwd(initial.dir)