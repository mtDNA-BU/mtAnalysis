#' mtReadVCF_mtc function
#'
#' Extract allele, freq and coverage datasets from the VCF file containing
#' results called by mitoCaller. The datasets are the input of mtAAF and
#' mtSummary functions.
#' @param VCF_file a string of path and name of VCF file to import
#' @param n_sample number of subjects
#' @return A list contains the three matrices of allele, freq and coverage
#' with dimensions of 16569 x N. Rows correspond to loci and columns correspond
#' to subjects. It contains subject ID as the column names.
#' @import data.table stringi
#' @export
#' @examples
#'


# VCF_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/Dissertation/VCF/OutputCHS.vcf"
mtReadVCF_mtc <- function(VCF_file, n_sample){

  dat_vcf <- fread(VCF_file, skip="CHROM", data.table=F)
  pos_var <- dat_vcf$POS
  pos_nonvar <- setdiff(c(1:.mtLength), dat_vcf$POS)
  allele_pos_nonvar <- .mtRef[pos_nonvar]
  allele_pos_nonvar_mat <-replicate(n_sample, allele_pos_nonvar)

  allele <- matrix("0", nrow=.mtLength, ncol=n_sample)
  freq <- matrix("1", nrow=.mtLength, ncol=n_sample)
  coverage <- matrix(NA_character_, nrow=.mtLength, ncol=n_sample)

  allele[pos_nonvar,] <- allele_pos_nonvar_mat


  # split the values of allele, freq and coverage from VCF
  dat_vcf_geno <- as.matrix(dat_vcf[,c((ncol(dat_vcf)-n_sample+1):ncol(dat_vcf))])
  mat_vcf_geno <- stri_split_fixed(dat_vcf_geno, ":", simplify = TRUE)

  x1k <- mat_vcf_geno[,1]
  x2k <- mat_vcf_geno[,2]
  x1k[x1k == "."] <- "0/0"
  x2k[x2k == "."] <- "0"


  x1k_m <- stri_split_fixed(x1k, "/", simplify = TRUE)
  class(x1k_m) <- "integer"
  x1k_m[is.na(x1k_m)] <- -1
  dta_vcf_alt_m <- stri_split_fixed(dat_vcf$ALT, pattern=",", simplify = TRUE)
  Both <- cbind(dat_vcf$REF, dta_vcf_alt_m )
  allele[pos_var, ] <- GetAlleles(Both, x1k_m, n_sample )
  rownames(allele) <- 1:.mtLength
  colnames(allele) <- names(dat_vcf)[c((ncol(dat_vcf)-n_sample+1):ncol(dat_vcf))]

  x2k_m <- stri_split_fixed(x2k, ",", simplify = TRUE)
  class(x2k_m) <- "numeric"
  x2k_m[is.na(x2k_m)] <- 0
  freq[pos_var, ] <- GetFreq(x2k_m)
  rownames(freq) <- 1:.mtLength
  colnames(freq) <- names(dat_vcf)[c((ncol(dat_vcf)-n_sample+1):ncol(dat_vcf))]


  coverage[pos_var, ] <- matrix(mat_vcf_geno[,5], ncol=n_sample)
  coverage[coverage == "."] <- NA_character_
  class(coverage) <- "integer"
  rownames(coverage) <- 1:.mtLength
  colnames(coverage) <- names(dat_vcf)[c((ncol(dat_vcf)-n_sample+1):ncol(dat_vcf))]
  #dim(coverage) <- c(nrow(Both), n_sample)


  out_list <- list("allele" = allele, "freq" = freq, "coverage" = coverage)

  out_list
}










