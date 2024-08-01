# qtlseq.R User Guide

## Table of contents
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Outputs](#outputs)



## Dependencies
### R libraries
- [tidyverse](https://www.tidyverse.org/)
- [argparse](https://cran.r-project.org/web/packages/argparse/index.html)
- [patchwork](https://patchwork.data-imaginist.com/)

#### Installation using conda

Create a QTL-seq specific environment with R. 
```
conda create -n qtlseq
conda activate qtlseq
```
Then install the packages using conda.
```
conda install -c conda-forge r-base r-essentials
conda install -c r r-tidyverse r-argparse r-patchwork r-locfit
```
## Usage

```
Rscript qtlseq.R -h

usage: qtlseq.R [-h] [-v <VCF>] [-b1 <CHR>] [-b2 <CHR>]
                   [-c <CHR> [<CHR> ...]] [-o <DIR>] [-w <INT>] [-ms <DOUBLE>]
                   [-rf <DOUBLE>] [-d <INT>] [-D <INT>] [--fig-width <DOUBLE>]
                   [--fig-height <DOUBLE>]

options:
  -h, --help            show this help message and exit
  -v <VCF>, --vcf <VCF>
                        VCF file. This VCF file must have AD field.
  -b1 <CHR>, --highbulk <CHR>
                        High bulk name.
  -b2 <CHR>, --lowbulk <CHR>
                        Low bulk name
  -o <DIR>, --out <DIR>
                        Output directory.                      
  -c <CHR> [<CHR> ...], --chrom <CHR> [<CHR> ...]
                        A list of chromosomes to be included in the analysis,
                        separated by space
  -w <INT>, --window <INT>
                        Window size (kb). [25000]
  -ms <DOUBLE>, --min-SNPindex <DOUBLE>
                        Cutoff of minimum SNP-index. [0.3]
  -rf <DOUBLE>, --ref-frq <DOUBLE>
                        Cutoff of reference allele frequency. Range will be
                        [rf] to [1 - rf]. [0.3]
  -d <INT>, --min-depth <INT>
                        Minimum depth of variants. [8]
  -D <INT>, --max-depth <INT>
                        Maximum depth of variants. [250]
  --fig-width <DOUBLE>  Width allocated in chromosome figure. [10]
  --fig-height <DOUBLE>
                        Height allocated in chromosome figure. [5]

```

### Parameters:
How reference allele frequency, SNP index for each bulk and delta SNP index are calculated:
![calc](https://github.com/Brycealong/QTL-analysis/blob/main/images/calc.png)

A note about window sizes:
The calculations are performed using the `locfit` function from the locfit package using a user defined window size and the degree of the polynomial set to zero. For a discussion about window size, we recommend reading Magwene et al. (2011). In general, **larger windows will produce smoother data**. The functions making these calculations are rather fast, so we recommend testing several window sizes for your data, and then deciding on the optimal size.

### Example:

```
Rscript qtlseq.R -v wheat-vcf/ALL.vcf.gz \
        -b1 Mutant123 \
        -b2 Wild123 \
        -o OUT_DIR \
        -c chr1 chr2 chr3 \
        -w 25000
```

### Possible Errors:

You may encounter `Out of vertex space` error when doing smoothing. Have the window size go higher will generally solve this issue. For more: https://github.com/bmansfeld/QTLseqr/issues/60

## Outputs
Files inside `OUT_DIR` are like below.
```
├── allchr.png
├── chr1.png
├── chr2.png
├── ...
├── SNPindex.tsv
├── SNPindex.filt.tsv
├── distribution.png
├── analysis.log
```
- check the results.
  + `SNPindex.tsv` : columns in this order.
    + **CHROM** - The chromosome this SNP is in 
    + **POS** - The position on the chromosome in nt 
    + **REF** - The reference allele at that position 
    + **ALT** - The alternate allele 
    + **DP.HIGH** - The read depth at that position in the high bulk 
    + **AD_REF.HIGH** - The allele depth of the reference allele in the high bulk 
    + **AD_ALT.HIGH** - The alternative allele depth in the high bulk  
    + **SNPindex.HIGH** - The calculated SNP-index for the high bulk 
    + Same as above for the low bulk 
    + **REF_FRQ** - The reference allele frequency as defined above 
    + **deltaSNP** - The $\Delta$(SNP-index) as defined above
  + `SNPindex.filt.tsv` : SNPs filtered with user-specified or default thresholds. One column `tricubeDeltaSNP` is added, which represents the smoothed deltaSNP values.
  + `allchr.png` : delta SNP index for all chromosomes
    - **dots** : variant
    - **<span style="color: red; ">RED line</span>** : smoothed delta SNP-index
    ![allchr](https://github.com/Brycealong/QTL-analysis/blob/main/images/allchr.png)
  + `chr1.png` : delta SNP index for one chromosome. Same for other chromosomes.
  ![6a](https://github.com/Brycealong/QTL-analysis/blob/main/images/6A.png)
  + `distribution.png`: distribution of reference allele frequency, read depths of each sample and SNP index of each sample. Adjust your thresholds using this graph.
  ![dis](https://github.com/Brycealong/QTL-analysis/blob/main/images/distribution.png)
  + `analysis.log`: log how many SNPs are filtered out on each parameter.

  
## Citation
- Hiroki Takagi, Akira Abe, Kentaro Yoshida, Shunichi Kosugi, Satoshi Natsume, Chikako Mitsuoka, Aiko Uemura, Hiroe Utsushi, Muluneh Tamiru, Shohei Takuno, Hideki Innan, Liliana M. Cano, Sophien Kamoun, Ryohei Terauchi (2013).  [QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations](https://doi.org/10.1111/tpj.12105). Plant journal 74:174-183.
