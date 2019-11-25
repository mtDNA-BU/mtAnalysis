#' mtdna test function
#'
#' Test function for the package
#'
#' @param head_file Head File
#' @param coverage_file Coverage File
#' @param allele_file Allele File
#' @param freq_file Frequency File
#'
#' @return character
#' @export
#' @examples
#' head_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/allele/allele.csv"
#' coverage_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/coverage/coverage.csv"
#' allele_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/allele/allele.csv"
#' freq_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/freq/freq.csv"
#' mtdna_test (head_file, coverage_file, allele_file, freq_file)
mtdna_test <- function( head_file, coverage_file, allele_file, freq_file){

  head <- scan(head_file, sep = ",", character(), nlines = 1)

}
