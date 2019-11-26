
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mtdnaANNO

<!-- badges: start -->

<!-- badges: end -->

Association analysis of heteroplasmic mtDNA mutations.

## Installation

You can install the development version of mtdnaANNO packag from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mtDNA-BU/mtdnaANNO")
```

## Example

This is an example how to run mtPAA() function

``` r
library(mtdnaANNO)

# Read input data
coverage_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/coverage/coverage.csv"
allele_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/allele/allele.csv"
freq_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/freq/freq.csv"

#head <- scan( file = allele_file, sep = ",", character(), nlines = 1, quiet = TRUE)
#coverage <- matrix( scan( file = coverage_file, sep = ",", character() ), ncol = length( head ), byrow = TRUE)
#allele <- matrix( scan( file = allele_file, sep = ",", character() ), ncol = length( head ), byrow = TRUE)
#freq <- matrix( scan( file = freq_file, sep = ",", character() ), ncol = length( head ), byrow = TRUE)

#PAA <- mtPAA (coverage, allele, freq)
#dim (PAA)
```
