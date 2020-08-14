
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ANNOmtDNA

### Association analysis of heteroplasmic mtDNA mutations.

**Authors:** Xianbang Sun (maintainer, <sxb3000@bu.edu>), Katia
Bulekova, Kaiyu Yan, Daniel Levy, Jun Ding, Jessica L. Fetterman, Chunyu
Liu<br> **Date:** “08/14/2020”

<!-- badges: start -->

<!-- badges: end -->

ANNOmtDNA is an R package for identifying and annotating mtDNA sequence
variations. This package allows to identify mtDNA sequence variations
with user-specified thresholds, and provide several plots to visualize
mtDNA sequence variations. It also annotates all mtDNA variations and
predicts functionality for mtDNA variants in the coding and tRNA
regions.

*Note:* if you use ANNOmtDNA in published research, please cite:

Xianbang Sun, Katia Bulekova, Kaiyu Yan, Daniel Levy, Jun Ding, Chunyu
Liu, Jessica L. Fetterman (2020) mtdnaANNO: an R package for
comprehensive annotation of mtDNA sequence variation

## Installation

You can install the development version of ANNOmtDNA package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("mtDNA-BU/ANNOmtDNA")
library(ANNOmtDNA)
```

## Standard workflow

Import the original allele, frequency and coverage datasets

``` r
input_path = "/path/to/input/directory/"
```

``` r
allele <- read.csv(paste0(input_path, "allele/allele_ANNOmtDNA.csv"), header = T, stringsAsFactors=FALSE,
                   colClasses = c("character"))
freq <- read.csv(paste0(input_path,"freq/freq_ANNOmtDNA.csv"), header=T, stringsAsFactors=FALSE,
                      colClasses = c("character"))
coverage <- read.csv(paste0(input_path,"coverage/coverage_ANNOmtDNA.csv"), header=T)
allele <- as.matrix(allele)
freq <- as.matrix(freq)
coverage <- data.matrix(coverage)
```

allele: a character matrix (16569 x N) provided by the user. Rows
correspond to loci and columns correspond to subjects. This matrix
contains N subjects with mtDNA sequencing data of 16569 loci. The matrix
must contain subject ID as the column names. “/” is used to delimited
different allele calls in a locus.

freq: a character matrix (16569 x N) provided by the user. Rows
correspond to loci and columns correspond to subjects. This matrix
contains the N subjects with mtDNA sequencing data of 16569 loci. The
matrix must contain subject ID as the column names. “/” is used to
delimited the allele fractions.

coverage: a numeric matrix (16569 x N). Rows correspond to loci and
columns correspond to subjects. This matrix contains the reads coverage
of the 16569 mtDNA loci for each subject. The matrix must contain the
subject ID as the column names.

#### Compute alternative allele fraction (AAF) by mtAAF function

``` r
AAF <- mtAAF(allele, freq)
```

Scatter plot of output of mtAAF function, each point is AAF of a subject
at a locus, the x axis is the mtDNA loci and the y axis is the range of
AAF (0-100%)

``` r
plot(AAF)
```

<img src="man/figures/README-scatterPlot-1.png" style="display: block; margin: auto;" />

Plot the mean coverage of each locus by plotCover function, the x axis
is the mtDNA loci and y axis is the mean coverage

``` r
plotCover(coverage)
```

<img src="man/figures/README-plotMeanCoverage-1.png" style="display: block; margin: auto;" />

Histogram of the mean coverage of each subject by histSampCov function

``` r
histSampCov(coverage)
```

<img src="man/figures/README-histMeanCoverage-1.png" style="display: block; margin: auto;" />

#### Summarize and annotate mtDNA mutations by mtSummary function

Specify the path to output the annotation file and histograms of mtDNA
mutation burden of subjects and number of mutations carried by each loci
for heteroplasmic and homoplasmic mutations

Users can specify lower bound and upper bound of the threshold by
thre.lower and thre.upper arguments. By default, thre.lower=0.03 and
thre.upper=0.97. That is, if \(0.03\leq AAF\leq 0.97\), it is a
heteroplasmic mutation; and if \(AAF>0.97\), it is a homoplasmic
mutation. Users can also choose different gene regions by loci argument.
For example, to choose tRNA region, set loci=“tRNA”. Users can also
choose types of mutations to be annotated by type argument. For example,
to choose heteroplasmic mutations, set type=“heter”, and choose
homoplasmic mutations, set type=“homo”. Run the mtSummary function, only
include loci of coding region, annotate for both of heteroplasmic and
homoplasmic mutations.

``` r
output_path <- "/output/dir/"
```

``` r
mtSum <- mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,
         path=output_path, type="both", study="ARIC")
```

Part of the output of annotated alleles

    #>     mtID ref_allele allele_var n_var n_heter n_homo mut_allele  Pos ref
    #> 766 3310          C        C/T     1       0      1          T 3310   2
    #>     Gene TypeMutation MissensMutation CodonPosition ProteinDomain
    #> 766   91            3            4636             1            10
    #>     dbSNP_150_id PolyPhen2 PolyPhen2_score SIFT SIFT_score CADD CADD_score
    #> 766           NA         1            0.11    2       0.34    1       2.45
    #>     CADD_phred_score
    #> 766            19.15

Output of histograms of summarized mutations

<img src="man/figures/README-mtHistograms-1.png" width="60%" height="60%" style="display: block; margin: auto;" /><img src="man/figures/README-mtHistograms-2.png" width="60%" height="60%" style="display: block; margin: auto;" /><img src="man/figures/README-mtHistograms-3.png" width="60%" height="60%" style="display: block; margin: auto;" /><img src="man/figures/README-mtHistograms-4.png" width="60%" height="60%" style="display: block; margin: auto;" />

Summary of the mean coverage of loci

``` r
mtSum$coverLoci 
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   145.9  5333.5  5572.2  5513.0  5727.4  6230.7
```

Summary of the mean coverage of subjects

``` r
mtSum$coverSubjects  
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   971.4  4334.7  5317.6  5513.0  6462.6 16415.2
```

Display loci of variation

``` r
mtSum$loci_var  
```

Summary of the heteroplasmic burden of subjects

``` r
mtSum$heter_burden_sum  
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.000   1.000   1.000   1.917   2.000 111.000
```

Summary of numbers of heteroplasmic mutations loci carried

``` r
mtSum$heter_loci_sum  
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#>    0.000    0.000    0.000    0.466    0.000 1877.000
```

Display loci of heteroplasmy

``` r
mtSum$loci_heter  
```

Total Number of heteroplasmic mutations

``` r
mtSum$heter_total  
#> [1] 7719
```

Summary of the homoplasmic burden of subjects

``` r
mtSum$homo_burden_sum  
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     0.0    13.0    27.0    25.1    34.0    96.0
```

Summary of numbers of homoplasmic mutations loci carried

``` r
mtSum$homo_loci_sum  
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#>    0.000    0.000    0.000    6.102    0.000 3975.000
```

Display loci of homoplasmy

``` r
mtSum$loci_homo  
```

Total Number of homoplasmic mutations

``` r
mtSum$homo_total  
#> [1] 101075
```

#### Annotate alternative alleles by mtAnno function

Generate alternative alleles to be annotated for mtDNA loci. It has two
columns: loci positions (“pos”) and “alleles” to be annotated.

``` r
anno <- as.data.frame(matrix(0, 10, 2))
colnames(anno)<- c("pos", "alleles")
anno$pos <- c(3311:3320)
anno$alleles <- c("A", "A", "T", "C", "A", "T", "G", "A", "C", "T")
```

Run the mtAnno function

``` r
mtAnno(anno=anno, path=output_path)
```

Part of the output of annotated alleles

    #>    pos alleles  Pos ref Gene  TypeMutation MissensMutation CodonPosition
    #> 1 3311       A 3311   C  ND1 Nonsynonymous             P2H             2
    #>            ProteinDomain dbSNP_150_id PolyPhen2 PolyPhen2_score    SIFT
    #> 1 Transmembrane; Helical           NA    benign            0.01 neutral
    #>   SIFT_score        CADD CADD_score CADD_phred_score
    #> 1       0.34 deleterious       2.45            19.12
