#' mtPAA function
#'
#' Calculate the Proportion of Alternative Allele of Each Locus for Each Individual
#'
#' @param coverage Coverage, numeric matrix
#' @param allele Allele dataset, numeric matrix
#' @param freq Frequency dataset, numeric matrix
#'
#' @return PAA, numeric matrix with dimensions of n subjects x 16569 SNPs
#' @export
#' @examples
#'
#' #Read input data
#' coverage_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/coverage/coverage.csv"
#' allele_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/allele/allele.csv"
#' freq_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/freq/freq.csv"
#'
#' head <- scan( file = allele_file, sep = ",", character(), nlines = 1, quiet = TRUE)
#' coverage <- matrix( scan( file = coverage_file, sep = ",", character() ), ncol = length( head ), byrow = TRUE)
#' allele <- matrix( scan( file = allele_file, sep = ",", character() ), ncol = length( head ), byrow = TRUE)
#' freq <- matrix( scan( file = freq_file, sep = ",", character() ), ncol = length( head ), byrow = TRUE)
#'
#' mtPAA (coverage, allele, freq)
mtPAA <- function( coverage, allele, freq ){


  if((sum(dim(coverage) != dim(allele) ) > 0) | (sum(dim(coverage) != dim(freq))>0)){
    stop("the coverage, allele and frequency should have same dimension")
  }

  # order the allele, freq and coverage datasets by ID
  allele <- allele[order(allele[,1]), ]
  freq <- freq[order(freq[,1]), ]
  coverage <- coverage[order(coverage[,1]), ]

  # record the IDs in the vector of subjectID, and remove the ID column from the original dataset
  subjectID <- allele[,1]
  allele <- allele[,-1]
  freq   <- freq[,-1]
  coverage<- coverage[,-1]

  # transpose the allele, freq and coverage datasets, so each column vector represents the sequence of a subject
  allele <- t(allele)
  freq <- t(freq)
  coverage <- t(coverage )


  ### to make the function be flexible for any loci, we need an additional parameter of loci position
  ### it is vector provided by users, with default of c(1:16569)
  loci <- c(1:16569)


  if (length(loci)!=(dim(allele)[1])){
    warning("the length of loci should be the same as the dataset")
  }

  if(length(loci)>16569){
    warning("loci should have no more than 16569")
  }

  AAF3.m <- array(0,rev(dim(allele)) )
  rownames(AAF3.m) <- subjectID

  for (m in 1:dim(allele)[2]) {

    m.allele   <- allele[,m]
    m.freq     <- freq[,m]

    complex.allele <- grep("/", m.allele)
    AAF2 <- array(NA, 16569)

    if(length(complex.allele)==0){
      AAF2 <- ifelse(m.allele == ref, 0, 1)
    }else{
      AAF2[-complex.allele] <- ifelse(m.allele[-complex.allele] == ref[-complex.allele], 0, 1)

      complex.all2 <- strsplit(m.allele[complex.allele],split="/")
      complex.freqsplit <- strsplit(m.freq[complex.allele],split="/")
      complex.freqsplit <- lapply(complex.freqsplit, as.numeric)
      complex.ref <- ref[complex.allele]
      AAF2[complex.allele] <- mapply(function(x, y, z){ pos <- which(x %in% z); if(length(pos) > 0) max(y[-pos]) else max(y)  } , x= complex.all2, y = complex.freqsplit, z = complex.ref)
    }

    AAF3.m[ m, ] <- AAF2
  }
  AAF3.m

}
