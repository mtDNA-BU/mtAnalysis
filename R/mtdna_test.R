#' mtdna test function
#'
#' Test function for the package
#'
#' @param coverage_file Coverage File
#' @param allele_file Allele File
#' @param freq_file Frequency File
#' @param path Path
#'
#' @return character
#' @export
#' @examples
#' coverage_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/coverage/coverage.csv"
#' allele_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/allele/allele.csv"
#' freq_file <- "/restricted/projectnb/mtdna-alcohol/Sun_Xianbang/ARIC/freq/freq.csv"
#' path<-"/rprojectnb2/mtdna-alcohol/Sun_Xianbang/Annotation/Output/"
#' #mtdna_test (coverage_file, allele_file, freq_file)
mtdna_test <- function( coverage_file, allele_file, freq_file, path){

  head  <- scan( file = allele_file, sep = ",", character(), nlines = 1, quiet = TRUE)
  coverage <- matrix( scan( file = coverage_file, sep = ",", character() ),
                      ncol = length( head ), byrow = TRUE)
  allele <- matrix( scan( file = allele_file, sep = ",", character() ),
                    ncol = length( head ), byrow = TRUE)
  freq <- matrix( scan( file = freq_file, sep = ",", character() ),
                  ncol = length( head ), byrow = TRUE)

  rownames(allele) <- allele[,1]
  rownames(freq)<-freq[,1]
  rownames(coverage)<-coverage[,1]

  # record the IDs in the vector of subjectID, and remove the ID column from the original dataset
  subjectID <- rownames(allele)
  allele <- allele[,-1]
  freq   <- freq[,-1]
  coverage<- coverage[,-1]

  #Transpose allele,  freq, and coverage file
  allele <- t(allele)
  freq <- t(freq)
  coverage <- t(coverage)
  class(coverage) <- "numeric"

  # order the allele, freq and coverage datasets by ID
  subjectID<-sort(subjectID)
  allele <- allele[ , subjectID]
  freq <- freq[ , subjectID]
  coverage <- coverage[ , subjectID]

  #================
  # mtAAF
  #================

  system.time(AAF <- mtAAF(allele, freq))

  #================
  # plot.mtDNAaaf
  #================

  system.time(plot(AAF))

  #================
  # plotCover
  #================

  plotCover(coverage)
  plotCover(coverage, "tRNA")

  #================
  # histSampCov
  #================

  histSampCov(coverage)
  histSampCov(coverage, "tRNA")

  #================
  # mtSummary
  #================
path = "/rprojectnb/mtdna-alcohol/Katia/"
  system.time(mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,loci="coding"
                        ,path=path, type="both"))

  system.time(mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,loci="coding"
                        ,path=path, type="both", study="ARIC"))

  system.time(mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,loci="coding"
                        ,path=path, type="heter", study="ARIC"))

  system.time(mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,loci="coding"
                        ,path=path, type="homo", study="ARIC"))


}
