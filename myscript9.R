#!usr/bin/env Rscript

# load the necessary libraries
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(argparse)
    library(patchwork)
  }
)

# Define command-line options

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--vcf",
                    help="VCF file. This VCF file must have AD field.",
                    metavar="<VCF>")
parser$add_argument("--highbulk",
                    help="high bulk name",
                    metavar="<CHR>")
parser$add_argument("--lowbulk", 
                    help="low bulk name",
                    metavar="<CHR>")
parser$add_argument("--chrom", type="character", nargs = "+",
                    help="a list of chromosomes to be included in the analysis, separated by space",
                    metavar="<CHR>")
parser$add_argument("-o", "--out", default = "out",
                    help="Output directory. Specified name can exist.",
                    metavar="<DIR>")
parser$add_argument("-w", "--window", default=2000, type="integer",
                    help="Window size (kb). [%(default)s]",
                    metavar="<INT>")
parser$add_argument("--min-SNPindex", default=0.3, type="double",
                    help="Cutoff of minimum SNP-index for clear results. [%(default)s]",
                    metavar="<DOUBLE>")
parser$add_argument("--ref-frq", default=0.3, type="double",
                    help="Cutoff of reference allele frequency. Range will be [%(default)s] to [1 - %(default)s].",
                    metavar="<DOUBLE>")
parser$add_argument("-d", "--min-depth", default=8, type="integer",
                    help="Minimum depth of variants which will be used. [%(default)s]",
                    metavar="<INT>")
parser$add_argument("-D", "--max-depth", default=250, type="integer",
                    help="Maximum depth of variants which will be used. [%(default)s]",
                    metavar="<INT>")
parser$add_argument("--fig-width", type = "double", default = 10.0, 
                    help = "Width allocated in chromosome figure. [%(default)s]",
                    metavar="<DOUBLE>")
parser$add_argument("--fig-height", type = "double", default = 5.0, 
                    help = "Height allocated in chromosome figure. [%(default)s]",
                    metavar="<DOUBLE>")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# Check for required arguments
if (is.null(args$vcf) || is.null(args$highbulk) || is.null(args$lowbulk)) {
  parser$print_help()
  stop("At least three arguments must be supplied (input VCF file , high bulk name, low bulk name).")
}

if (is.null(args$chrom))
  message("Chromosome list not provided. Using all the chromosomes in VCF file...\n",
  "Highly recommend using `bcftools query -f'%CHROM\\n' [file.vcf.gz] | uniq -c` to check any chromosomes unwanted.")

importFromVCF <- function(file,
                          highBulk,
                          lowBulk,
                          chromList = NULL){
  extractAD <- function(format, sample){
    # split the format and sample strings
    format.parts <- unlist(strsplit(format, ":"))
    sample.parts <- unlist(strsplit(sample, ":"))
    
    # find the position of 'AD' in the format
    ad.idx <- base::which(format.parts == "AD")
    # find the 'AD' value from the sample
    ad.val <- sample.parts[ad.idx]
    return(ad.val)
  }
  
  message("Reading and processing VCF file...")
  SNPset <- readr::read_delim(file, comment = "##", show_col_types = F) %>%
    dplyr::rename(CHROM = "#CHROM")
  
  if (!all(c(highBulk, lowBulk) %in% colnames(SNPset))){
    # whichnot <-
    #   paste(chromList[base::which(!chromList %in% unique(SNPset$CHROM))], collapse = ', ')
    stop("Bulk names incorrect. Use `bcftools query -l [file.vcf.gz]` to check the bulk names.")
  } else {
    #Keep only wanted chromosomes
    if (!is.null(chromList)) {
      if (!all(chromList %in% unique(SNPset$CHROM))) {
        whichnot <-
          paste(chromList[base::which(!chromList %in% unique(SNPset$CHROM))], collapse = ', ')
        stop(paste0("The following are not true chromosome names: ", whichnot))
      } else {
        message("Removing the following chromosomes: ", paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in% chromList], collapse = ", "))
        SNPset <- SNPset[SNPset$CHROM %in% chromList, ]
      }
    } else {
      SNPset <- SNPset %>%
        dplyr::select(CHROM, POS, REF, ALT, FORMAT, all_of(c(highBulk, lowBulk))) %>%
        dplyr::filter(!grepl(",", ALT)) %>% # remove multi-allelic
        dplyr::mutate(AD.HIGH = purrr::map2_chr(FORMAT, .data[[highBulk]], extractAD)) %>%
        dplyr::mutate(AD.LOW = purrr::map2_chr(FORMAT, .data[[lowBulk]], extractAD)) %>%
        dplyr::select(!FORMAT) %>%
        dplyr::select(!all_of(c(highBulk, lowBulk))) %>%
        tidyr::separate(
          col = AD.HIGH,
          into = c("AD_REF.HIGH", "AD_ALT.HIGH"),
          sep = ",",
          convert = T
        ) %>%
        tidyr::separate(
          col = AD.LOW,
          into = c("AD_REF.LOW", "AD_ALT.LOW"),
          sep = ",",
          convert = T
        ) %>%
        dplyr::mutate(
          DP.HIGH = AD_REF.HIGH + AD_ALT.HIGH,
          DP.LOW = AD_REF.LOW + AD_ALT.LOW
        ) %>%
        dplyr::mutate(
          SNPindex.HIGH = AD_ALT.HIGH / DP.HIGH,
          SNPindex.LOW = AD_ALT.LOW / DP.LOW,
          REF_FRQ = (AD_REF.HIGH + AD_REF.LOW) / (DP.HIGH + DP.LOW),
          deltaSNP = SNPindex.HIGH - SNPindex.LOW
        ) %>%
        dplyr::mutate_all(~ ifelse(is.na(.), 0, .))
    }
    SNPset
  }
}

filterSNPs <- function(SNPset,
                       minSNPindex,
                       refAlleleFreq,
                       minSampleDepth,
                       maxSampleDepth,
                       verbose = TRUE) {
  
  org_count <- nrow(SNPset)
  count <- nrow(SNPset)
  
  # Filter by total reference allele frequency
  if (!is.null(minSNPindex)) {
    if (verbose) {
      message(
        "Filtering by either sample SNP index >= ",
        minSNPindex
      )
    }
    SNPset <- SNPset %>% 
      dplyr::filter(SNPindex.HIGH >= minSNPindex | SNPindex.LOW >= minSNPindex)
    if (verbose) {
      message("...Filtered ", count - nrow(SNPset), " SNPs")
    }
    count <- nrow(SNPset)
  }

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
    SNPset <- dplyr::filter(SNPset, SNPset$REF_FRQ <= 1 - refAlleleFreq &
                              SNPset$REF_FRQ >= refAlleleFreq)
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
  
  # Filter by max read depth in either sample
  if (!is.null(maxSampleDepth)) {
    if (verbose) {
      message("Filtering by per sample read depth: DP <= ",
              maxSampleDepth)
    }
    SNPset <-
      dplyr::filter(SNPset, DP.HIGH <= maxSampleDepth &
                      DP.LOW <= maxSampleDepth)
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



tricubeStat <- function(POS, Stat, windowSize)
{
  if (windowSize <= 0)
    stop("A positive smoothing window is required")
  stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0)), POS)
}


runQTLseqAnalysis <- function(SNPset, windowSize){
  
  message("Calculating tricube smoothed delta SNP index...")
  return(SNPset %>%
           dplyr::group_by(CHROM) %>%
           dplyr::mutate(tricubeDeltaSNP = tricubeStat(POS = POS, Stat = deltaSNP, windowSize))
  )
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


plotQTLManh <-
  function(SNPset){
    
    data_cum <- SNPset |>
      group_by(CHROM) |>
      summarise(max.POS = max(POS)) |>
      mutate(POS.add = lag(cumsum(max.POS), default = 0)) |>
      select(CHROM, POS.add)
    
    SNPset <- SNPset |>
      inner_join(data_cum, by = "CHROM") |>
      mutate(POS.cum = POS + POS.add)
    
    axis_set <- SNPset |>
      group_by(CHROM) |>
      summarize(center = mean(POS.cum))
    
    p <- ggplot(SNPset) +
      geom_point(aes(x = POS.cum, y = deltaSNP,
                     color = forcats::as_factor(CHROM)), alpha = 0.75, size = 0.25) +
      geom_line(aes(x = POS.cum, y = tricubeDeltaSNP), color = "red", linewidth = 1) +
      scale_x_continuous(
        label = axis_set$CHROM,
        breaks = axis_set$center
      ) +
      scale_color_manual(values = rep(
        c("blue4", "green4"),
        unique(length(axis_set$CHROM))
      )) + 
      ggplot2::ylim(-1,1) +
      labs(
        x = NULL,
        y = expression(Delta * " (SNP-index)")
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, size = 8, vjust = 0.5),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
      )
    p
  }

plotQTLStats <-
  function(SNPset, chr) {
    
    chrset <- SNPset %>% dplyr::filter(CHROM == chr)
    p <- ggplot2::ggplot(data = chrset) +
      ggplot2::geom_point(ggplot2::aes(POS, deltaSNP), size = 0.5, color = "blue") +
      ggplot2::geom_line(ggplot2::aes(POS, tricubeDeltaSNP), linewidth = 1, color = "red") +
      ggplot2::ylim(-1,1) +
      ggplot2::scale_x_continuous(breaks = seq(from = 0, to = max(chrset$POS), by = 10^(floor(log10(max(chrset$POS))))), 
                                  labels = format_genomic(), 
                                  name = paste0(chr, " (", format_prefix(max(chrset$POS)), ")")) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = element_blank()) +
      ggplot2::ylab(expression(Delta * " (SNP-index)"))
      
  }
message("Creating output directory ", args$out, "...")
dir.create(args$out, recursive = T)
df <- importFromVCF(
  args$vcf, 
  highBulk = args$highbulk, 
  lowBulk = args$lowbulk, 
  chromList = args$chrom
)
p1 <- ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))
p2 <- ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.HIGH)) +
  xlim(0,1)
p3 <- ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.LOW)) +
  xlim(0,1)
p4 <- ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH)) +
  xlim(0,1000)
p5 <- ggplot(data = df) +
  geom_histogram(aes(x = DP.LOW)) +
  xlim(0,1000)

p <- p1/(p2|p3)/(p4|p5)

ggplot2::ggsave(file.path(args$out, "distribution.png"), p, 
                width = args$fig_width, height = args$fig_height)

readr::write_tsv(df, file = file.path(args$out, "SNPindex.tsv"))
df_filt <- filterSNPs(
  df,
  minSNPindex = args$min_SNPindex,
  refAlleleFreq = args$ref_frq,
  minSampleDepth = args$min_depth,
  maxSampleDepth = args$max_depth
)
df_filt <- runQTLseqAnalysis(df_filt, windowSize = args$window * 1000)
readr::write_tsv(df_filt, file = file.path(args$out, "SNPindex.filt.tsv"))
ggsave(file.path(args$out, "allchr.png"), plotQTLManh(df_filt), width = args$fig_width, height = args$fig_height)
for (chr in unique(df_filt$CHROM)) {
  p <- plotQTLStats(SNPset = df_filt, chr = chr)
  ggsave(file.path(args$out, paste0(chr, ".png")), p, width = args$fig_width, height = args$fig_height)
}
