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

  head <- scan( file = head_file, sep = ",", character(), nlines = 1)
  coverage <- matrix( scan( file = coverage_file, sep = ",", character() ),
                      ncol = length( head ), byrow = TRUE)
  allele <- matrix( scan( file = allele_file, sep = ",", character() ),
                    ncol = length( head ), byrow = TRUE)
  freq <- matrix( scan( file = freq_file, sep = ",", character() ),
                  ncol = length( head ), byrow = TRUE)


  if((sum(dim(coverage) != dim(allele) ) > 0) | (sum(dim(coverage) != dim(freq))>0)){
    stop("the coverage, allele and frequency should have same dimension")
  }

  # order the allele, freq and coverage datasets by ID
  allele <- allele[order(allele[,1]), ]
  freq <- freq[order(freq[,1]), ]
  coverage <- coverage[order(coverage[,1]), ]

  # record the IDs in the vector of subjectID, and remove the ID column from the original dataset
  subjectID <- allele[,1]
  allele2 <- allele[,-1]
  freq2   <- freq[,-1]
  coverage2<- coverage[,-1]

  ### to make the function be flexible for any loci, we need an additional parameter of loci position
  ### it is vector provided by users, with default of c(1:16569)
  loci <- c(1:16569)

  if (length(loci)!=(dim(allele2)[2])){
    warning("the length of loci should be the same as the dataset")
  }

  if(length(loci)>16569){
    warning("loci should have no more than 16569")
  }




}
