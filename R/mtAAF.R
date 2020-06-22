#' mtAAF function
#'
#' Derivation of alternative allele fraction
#'
#' @param allele data frame (N x 16569) provided by the user. This data frame contains
#' N subjects with mtDNA sequencing data of 16569 loci. The data frame contains subject ID as
#' the row names, and the allele calls of all mtDNA loci for each subject as the columns.
#' “/” is used to delimited different allele calls in a locus.
#' @param freq data frame (N x 16569) provided by the user. This data frame contains
#' N subjects with mtDNA sequencing data of 16569 loci. The data frame contains subject ID as
#' the row names, and the allele fractions of the called alleles for each subject as the columns.
#' “/” is used to delimited the allele fractions.
#' @return AAF, numeric matrix (N x 16569). It contains subject ID as the row names,
#' and the AAF of all 16569 mtDNA loci for each subject. It will also generate a scatter plot
#' of AAF based on user’s choice.
#' @import graphics
#' @export
#' @examples
#'
#' #Read input data
#' allele_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/allele/allele.csv"
#' freq_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/freq/freq.csv"
#'
#' #head <- scan( file = allele_file, sep = ",", character(), nlines = 1, quiet = TRUE)
#' #allele   <- matrix( scan( file = allele_file,
#' #                          sep = ",", character() ),
#' #                    ncol = length( head ), byrow = TRUE)
#' #freq     <- matrix( scan( file = freq_file,
#' #                          sep = ",", character() ),
#' #                    ncol = length( head ), byrow = TRUE)
#'
#' #mtPAA ( allele, freq)
mtAAF <- function( allele, freq ){

  if(!is.character(allele) | !is.character(freq)){
    stop("the mode of allele and freq must be character")
  }

  if(dim(allele)[1]!=16569 | dim(freq)[1]!=16569){
    stop("the allele and frequency should have 16569 rows")
  }

  if((sum(dim(allele) != dim(freq) ) > 0)){
    stop("the allele and frequency should have the same dimension")
  }

  # record the IDs in the vector of subjectID
  subjectID <- colnames(allele)

  # order the allele, freq datasets by ID
  subjectID<-sort(subjectID)
  allele <- allele[ , subjectID ]
  freq <- freq[ , subjectID]

  AAF3.m <- array(0, dim(allele) )
  complex.allele <- which(allele != .mtRef)
  AAF3.m [complex.allele] <- 1

  complex.allele <- complex.allele[ grepl( "/",  allele[complex.allele])]
  complex.all2 <- strsplit( allele[complex.allele], split="/" )

  complex.freqsplit <- strsplit( freq[complex.allele], split="/")
  complex.freqsplit <- lapply( complex.freqsplit, as.numeric)

  complex.ref <- .mtRef[ (complex.allele -1) %% length(.mtRef) + 1]

  AAF3.m[complex.allele] <- mapply(function(x, y, z){
    pos <- which(x == z);
    if(length(pos) > 0) 1-y[pos] else 1
  } ,
  x= complex.all2,
  y = complex.freqsplit,
  z = complex.ref)

  AAF3.m[is.na(AAF3.m)]<-0


  colnames(AAF3.m) <- subjectID
  rownames(AAF3.m) <- c(1:16569)
  class(AAF3.m) <- c('matrix', 'mtDNAaaf')
  AAF3.m

}
