# store the current directory
initial.dir <- getwd()
# change to the new directory
setwd("/home/bryce/Desktop/QTL_analysis/")
# # load the necessary libraries
library(dplyr)
library(tidyr)
# library(ggplot2)
# library(patchwork)
library(vcfR)

#load the package
library(QTLseqr)


#Set sample and file names
HighBulk <- "N_CD_188"
LowBulk <- "N_CH_188"
file <- "demo2/demo.2.vcf"

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- expand.grid(ChrNum = 1:7, ChrLetter = c("A", "B", "D")) %>% 
  mutate(Chrom = paste0("chr", ChrNum, ChrLetter)) %>% 
  pull(Chrom)

importFromVCF <- function(file,
                          highBulk = character(),
                          lowBulk = character(),
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
  
  SNPset <- tidy_gt %>%
    dplyr::filter(Indiv == LowBulk) %>% dplyr::select(-Indiv) %>%
    dplyr::left_join(dplyr::select(dplyr::filter(tidy_gt, Indiv == HighBulk),-Indiv),
                     by = "Key",
                     suffix = c(".LOW", ".HIGH")) %>%
    tidyr::separate(
      col = "AD.LOW",
      into = c("AD_REF.LOW", "AD_ALT.LOW"),
      sep = ",",
      # extra = "merge",
      convert = TRUE
    ) %>%
    tidyr::separate(
      col = "AD.HIGH",
      into = c("AD_REF.HIGH", "AD_ALT.HIGH"),
      sep = ",",
      # extra = "merge", 
      convert = TRUE
    ) %>%
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
  as.data.frame(SNPset)
}

#Import SNP data from file
df <-
  importFromVCF(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )

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
                      SNPset$DP.LOW >= minSampleDepth)
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

#Filter SNPs based on some criteria
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
    minGQ = 99
  )

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
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

#export summary CSV
getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")