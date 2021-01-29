#' mtAAF function
#'
#' Derivation of alternative allele fraction (AAF)
#'
#' @param allele a character matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains N subjects with mtDNA sequencing data of 16569 loci.
#' The matrix must contain subject ID as the column names.
#' "/" is used to delimited different allele calls in a locus.
#' @param freq a character matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains the N subjects with mtDNA sequencing data of 16569
#' loci. The matrix must contain subject ID as the column names.
#' "/" is used to delimited the allele fractions.
#' @param method the method to compute AAF for the case of multiple alternative
#' alleles. The default method "maxAA" computes AAF as the maximum of
#' frequencies of corresponding alternative alleles; and the alternative method
#' "allAA" computes AAF as 1 minus the frequency of reference allele
#' @return AAF, a numeric matrix (16569 x N).
#' Rows correspond to loci and columns correspond to subjects.
#' It contains subject ID as the column names, and the AAFs of all 16569 mtDNA
#' loci for each subject.
#' @import graphics
#' @export
#' @examples
#'
#'\dontrun{
#' ## Read input data
#' allele_file <- "allele.csv"
#' freq_file   <- "freq.csv"
#'
#' allele <- as.matrix( read.csv(file=allele_file, sep=",", header=FALSE))
#' freq   <- as.matrix( read.csv(file=freq_file, sep=",", header=FALSE))
#'
#' aaf=mtAAF(allele, freq)
#'}
#'
mtAAF <- function(allele, freq, method="maxAA"){

    if( !is.character(allele) | !is.character(freq) ){
        stop("the mode of allele and freq must be character")
    }

    if(dim(allele)[1] != .mtLength | dim(freq)[1] != .mtLength){
        stop("the allele and frequency should have 16569 loci (rows)")
    }

    if((sum(dim(allele) != dim(freq)) > 0)){
        stop("the allele and frequency should have the same dimension")
    }

    subjectID <- colnames(allele)
    subjectID<-sort(subjectID)
    allele <- allele[ , subjectID]
    freq <- freq[ , subjectID]

    AAF3.m <- array(0, dim(allele) )
    complex.allele <- which(allele != .mtRef)
    AAF3.m[complex.allele] <- 1

    complex.allele <- complex.allele[ grepl( "/", allele[complex.allele])]
    complex.all2 <- strsplit( allele[complex.allele], split="/" )

    complex.freqsplit <- strsplit( freq[complex.allele], split="/")
    complex.freqsplit <- lapply( complex.freqsplit, as.numeric)

    complex.ref <- .mtRef[ (complex.allele -1) %% length(.mtRef) + 1]

    if(method=="allAA"){
    AAF3.m[complex.allele] <- mapply(function(x, y, z){
        pos <- which(x == z);
        if(length(pos) > 0) 1-y[pos] else 1
    },
    x=complex.all2,
    y=complex.freqsplit,
    z=complex.ref)
    }else if(method=="maxAA"){
        AAF3.m[complex.allele] <- mapply(function(x, y, z){
            pos <- which(x == z);
            if(length(pos) > 0) max(y[-pos]) else max(y)
        },
        x=complex.all2,
        y=complex.freqsplit,
        z=complex.ref)
    }

    AAF3.m[is.na(AAF3.m)] <- 0
    AAF3.m[AAF3.m>1] <- 1
    AAF3.m[AAF3.m<0] <- 0

    colnames(AAF3.m) <- subjectID
    rownames(AAF3.m) <- seq_len(.mtLength)
    class(AAF3.m) <- c('matrix', 'mtDNAaaf')
    AAF3.m

}
