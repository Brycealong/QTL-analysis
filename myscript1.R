# store the current directory
initial.dir <- getwd()
# change to the new directory
setwd("/home/bryce/Desktop/QTL_analysis/")
# load the necessary libraries
library(dplyr)
library(ggplot2)
library(patchwork)
# load the dataset
vcf <- read.delim("vcf/4sample.all.parsed.vcf")
# look at the chromosomes
vcf %>% count(X.CHROM)
# remove chromosome 'Un'
vcf <- vcf %>% filter(X.CHROM != "Un")
# look at the ref allele summary
vcf %>% count(REF, sort = TRUE)
# look at the alt allele summary
vcf %>% count(ALT, sort = TRUE)
# remove SNPs that's not bi-allele
vcf_bi <- vcf %>% filter(!grepl(",", vcf$ALT))
# remove SNPs without 'AD' information
vcf_bi_ad <- vcf_bi %>% filter(grepl("AD", vcf_bi$FORMAT))

## Calculate the SNP index using the 'AD' information -----------
# make a copy of the filtered dataframe
df <- vcf_bi_ad
# function to extract the value from a sample column
calc_snpidx <- function(format, sample){
  # split the format and sample strings
  format.parts <- unlist(strsplit(format, ":"))
  sample.parts <- unlist(strsplit(sample, ":"))
  
  # find the position of 'AD' in the format
  ad.idx <- which(format.parts == "AD")
  # find the 'AD' value from the sample
  ad.val <- sample.parts[ad.idx]
  # split the value by comma and convert to numeric
  ad.nums <- as.numeric(unlist(strsplit(ad.val, "[,]")))
  
  # calculate and return the SNP index
  return(ad.nums[2] / sum(ad.nums))
}

# identify the sample columns
sample.cols <- names(df)[(which(names(df) == "FORMAT") + 1):ncol(df)]

# apply the function to each sample column and store the results in a new column
for (sample.col in sample.cols){
  df[[paste0(sample.col, ".snpidx")]] <- 
    mapply(calc_snpidx, df$FORMAT, df[[sample.col]])
}

# print(df)

# save the dataframe with indices
write.table(df, file = "KN1.T1-1.581.indices.vcf", sep = "\t", quote = FALSE)

# inspect the distribution of snp indices
df <- read.delim('KN1.T1-1.581.indices.vcf')
summary(df)
p1 <- ggplot(data = df, aes(x = M_22W581.2.snpidx)) + 
  stat_density()

p2 <- ggplot(data = df, aes(x = W_22W581.snpidx)) + 
  stat_density()

ggsave("idx_distribution.png", p1/p2, width = 5, height = 5)

# the correlation of two indices
# ggplot(df, aes(x = M_22W581.2.snpidx, y = W_22W581.snpidx)) +
#   geom_point()

## Filter SNPs ---------
# remove NAs
df <- na.omit(df)
# remove SNPs with SNP indices < 0.3 for both samples 
df <- df %>% filter(!(M_22W581.2.snpidx < 0.3 & W_22W581.snpidx < 0.3))

write.table(df, file = "KN1.T1-1.581.indices.filtered.vcf", sep = "\t", quote = FALSE)

## Sliding window analysis -----------
sliding_window_avg <- 
  function(df, window_size = 4e+07, increment_size=2e+05, 
           chr, idx_col){
    df <- df %>% filter(X.CHROM == chr)
    a <- min(df$POS)
    b <- a + window_size
    
    results <- data.frame(midpoint = numeric(), avg_val = numeric())
    
    # initialize start and end indices for the sliding window
    start_idx <- 1
    end_idx <- 1
    
    # iterate through the dataframe using the sliding window
    while (b <= max(df$POS)){
      # move the start_idx
      while (df$POS[start_idx] < a){
        start_idx <- start_idx + 1
      }
      # move the end_idx to include all elements up to 'b'
      while (df$POS[end_idx] <= b){
        end_idx <- end_idx + 1
      }
      # filter the dataframe using the sliding window 
      df1 <- df[start_idx:(end_idx - 1), ]
      
      if (nrow(df1) >= 10){
        val <- mean(df1[[idx_col]])
        results <- rbind(results, data.frame(midpoint = (a + b) / 2, avg_val = val))
      }
      # increment
      a <- a + increment_size
      b <- a + window_size
    }
    return(results)
  }

## Plot SNP index against position for each chromosome --------
# retrieve dataframe from file generated at previous step
df <- read.delim(file = 'KN1.T1-1.581.indices.filtered.vcf')
p1 <- ggplot(data = df %>% filter(X.CHROM == "1B"), 
             aes(x = POS/1e+06, y = M_22W581.2.snpidx)) +
  geom_point(color = "blue") +
  geom_hline(aes(yintercept = 0.5), linetype = 2) + 
  geom_line(data = sliding_window_avg(df, chr = "1B", idx_col = "M_22W581.2.snpidx"), 
            aes(x = midpoint/1e+06, y = avg_val), 
            color = "red") +
  labs(x = "Chromosome position (Mb)", y = "SNP index",
       title = "Chromosome 1B"
       )
# ggsave("plot0717_1.png", width = 5, height = 3)

p2 <- ggplot(data = df %>% filter(X.CHROM == "1B"), 
       aes(x = POS/1e+06, y = W_22W581.snpidx)) +
  geom_point(color = "blue") +
  geom_hline(aes(yintercept = 0.5), linetype = 2) + 
  geom_line(data = sliding_window_avg(df, chr = "1B", idx_col = "W_22W581.snpidx"), 
            aes(x = midpoint/1e+06, y = avg_val), 
            color = "red") +
  labs(x = "Chromosome position (Mb)", y = "SNP index",
       title = "Chromosome 1B"
       )
# ggsave("plot0717_1.png", width = 5, height = 3)
p <- p1 / p2
ggsave("chr1b.png", p, width = 5, height = 5)
# Manhattan plot ---------------
manhattan.plot <- function(chr, pos, pvalue, 
                           sig.level=NA, annotate=NULL, ann.default=list(),
                           should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                           xlab="Chromosome", ylab=expression(-log[10](p-value)),
                           col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()
  
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}

