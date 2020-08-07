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
#'
#'\dontrun{
#' #Read input data
#' coverage_file <- "coverage.csv"
#' allele_file <- "allele.csv"
#' freq_file   <- "freq.csv"
#'
#' coverage <- as.matrix( read.csv(file=coverage_file, sep=",", header=FALSE))
#' allele <- as.matrix( read.csv(file=allele_file, sep=",", header=FALSE))
#' freq   <- as.matrix( read.csv(file=freq_file, sep=",", header=FALSE))
#'
#' mtdna_test(coverage_file=coverage, allele_file=allele, freq_file=freq)
#'}
#'
mtdna_test <- function( coverage_file, allele_file, freq_file, path){

    head  <- scan( file=allele_file, sep=",", character(), nlines=1, quiet=TRUE)
    coverage <- matrix( scan( file=coverage_file, sep=",", character() ),
                        ncol=length( head ), byrow=TRUE)
    allele <- matrix( scan( file=allele_file, sep=",", character() ),
                      ncol=length( head ), byrow=TRUE)
    freq <- matrix( scan( file=freq_file, sep=",", character() ),
                    ncol=length( head ), byrow=TRUE)

    rownames(allele) <- allele[, 1]
    rownames(freq) <- freq[, 1]
    rownames(coverage) <- coverage[, 1]

    ## Record the IDs in the vector of subjectID, and remove the ID column from
    ## the original dataset
    subjectID <- rownames(allele)
    allele <- allele[, -1]
    freq   <- freq[, -1]
    coverage<- coverage[, -1]

    ## Transpose allele,  freq, and coverage file
    allele <- t(allele)
    freq <- t(freq)
    coverage <- t(coverage)
    class(coverage) <- "numeric"

    ## Order the allele, freq and coverage datasets by ID
    subjectID<-sort(subjectID)
    allele <- allele[, subjectID]
    freq <- freq[, subjectID]
    coverage <- coverage[, subjectID]

    ## mtAAF

    system.time(AAF <- mtAAF(allele, freq))

    ## plot.mtDNAaaf

    system.time(plot(AAF))

    ## plotCover

    plotCover(coverage)
    plotCover(coverage, "tRNA")

    ## histSampCov

    histSampCov(coverage)
    histSampCov(coverage, "tRNA")

    ## mtSummary

    path <- "/rprojectnb/mtdna-alcohol/Katia/"

    system.time(mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,
                          loci="coding", path=path, type="both", study="ARIC"))
    system.time(mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,
                          loci="coding", path=path, type="heter", study="ARIC"))
    system.time(mtSummary(aaf=AAF, allele=allele, freq=freq, coverage=coverage,
                          loci="coding", path=path, type="homo", study="ARIC"))


}
