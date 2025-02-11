#!usr/bin/env Rscript
Args <- commandArgs(trailingOnly = TRUE)

if (length(Args) < 3){
  stop("At least three argument must be supplied (input vcf) (high bulk name) (low bulk name)", call. = FALSE)
} else if (length(Args) == 3){
  # default output file
  message("using the default output file name 'out.png' ...")
  Args[4] <- "out.png"
}
# # store the current directory
# initial.dir <- getwd()
# # change to the new directory
# setwd("/home/bryce/Desktop/QTL_analysis/")
# load the necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
# library(patchwork)
library(vcfR)

#load the package
#library(QTLseqr)


#Set sample and file names
file <- Args[1]
HighBulk <- Args[2]
LowBulk <- Args[3]

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- expand.grid(ChrNum = 1:7, ChrLetter = c("A", "B", "D")) %>% 
  mutate(Chrom = paste0("chr", ChrNum, ChrLetter)) %>% 
  pull(Chrom)

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
                          chromList = NULL) {
  
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
  if (!is.null(chromList)) {
    message("Removing the following chromosomes: ", paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in% chromList], collapse = ", "))
    SNPset <- SNPset[SNPset$CHROM %in% chromList, ]
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
  if (!missing(refAlleleFreq)) {
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
  
  if (!missing(filterAroundMedianDepth)) {
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
  
  if (!missing(minTotalDepth)) {
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
  
  if (!missing(maxTotalDepth)) {
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
  if (!missing(minSampleDepth)) {
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
  if (!missing(minGQ)) {
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
  if (!missing(depthDifference)) {
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
  
  # #Filter SNP Clusters
  # if (!is.null(SNPsInCluster) & !is.null(ClusterWin)) {
  #     tmp <- which(diff(SNPset$POS, SNPsInCluster-1) < ClusterWin)
  # message("...Filtered ", count - nrow(SNPset), " SNPs")
  # count <- nrow(SNPset)
  # }
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
           var = "nSNPs",
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

#Set sample and file names
HighBulk <- "WT8214DPool"
LowBulk <- "MT8214UPool"
file <- "TC-CQZ-20220810-001-LJ008-pepper-vcf/ALL.vcf.gz"

HighBulk <- c("HH", "HH-1", "HH-2", "HH-3")
LowBulk <- c("LH", "LH-1", "LH-2", "LH-3", "LH_6")
file <- "demo1/demo.vcf.gz"

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- expand.grid(ChrNum = 1:7, ChrLetter = c("A", "B", "D")) %>% 
  mutate(Chrom = paste0("chr", ChrNum, ChrLetter)) %>% 
  pull(Chrom)

#Import SNP data from file
df <-
  importFromVCF(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )

#Filter SNPs based on some criteria
df_filt <-
  filterSNPs(
    SNPset = df,
    maxTotalDepth = 250,
    minSampleDepth = 7,
    #minGQ = 99
  )

ggplot(data = df) + 
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
  xlim(0,1000)

ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
  SNPset = df_filt,
  windowSize = 4e7,
  popStruc = "F2",
  bulkSize = c(25, 25),
  replications = 10000,
  intervals = c(95, 99)
)

#Plot
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

ggplot2::ggsave(Args[4], width = 8, height = 5)
#export summary CSV
# getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")

# change back to the original directory
# setwd(initial.dir)