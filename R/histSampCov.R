#' histSampCov function
#'
#' Produce a Histogram of Mean Coverage Across mtDNA Loci for Subjects
#'
#' @param coverage a numeric matrix (16569 x N).
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains the reads coverage of the 16569 mtDNA loci for each
#' subject. The matrix must contain the subject ID as the column names.
#' @param loci one of the following to specify mtDNA loci:
#' 1. a numeric vector (default is c(1:16569)) of mitochondrial DNA loci,
#' 2. a character string for the regions(e.g. "coding", "tRNA", "Dloop", …).
#' @param type a string: "median" to calculate median coverage at individual
#' level, "mean" to calculate mean coverage at individual level
#' @param ... arguments to be passed to plot function.
#' @export
#' @examples
#'
#'
#'\dontrun{
#' ## Read input data
#' coverage_file <- "coverage.csv"
#' coverage <- as.matrix(read.csv(file=coverage_file, sep=",", header=FALSE))
#'
#' histSampCov(coverage, loci="coding")
#'
#'}
#'
histSampCov <- function(coverage, loci=seq_len(.mtLength), type="median", ...) {

    if(length(type)!=1){
        stop("type must be a string of median or mean")
    }else if(type!="median" & type!="mean"){
        stop("type must be a string of median or mean")
    }

    if(dim(coverage)[1] != .mtLength){
        stop("the coverage should have 16569 loci (rows)")
    }

    if(is.numeric(loci) & !all(loci %in% seq_len(.mtLength))){
        stop("loci should be a subset of 1:16569")
    }

    ## assign loci value when the specified loci are strings (e.g. coding)
    if (is.character(loci)){
        if( loci == "coding"){
            loci <- .loci.coding
        }else if(loci == "tRNA"){
            loci <- .loci.tRNA
        }else if( loci == "RNR1"){
            loci <- .loci.RNR1
        }else if( loci == "RNR2"){
            loci <- .loci.RNR2
        }else {
            stop("loci name must be one of: coding, tRNA, RNR1, RNR2")
        }
    }

    coverage <- coverage[loci, ]
    rownames(coverage) <- as.character(loci)


    if(type=="median"){

        ## histogram of median coverage of subjects across loci
        cov_sub <- apply(coverage, 2, FUN = median)
        cov_sub_hist <- as.data.frame(cov_sub)

        h <- hist(cov_sub_hist[, 1],
                  xlab="Median coverage of subjects",
                  main="Histogram of Median Coverage of Subjects", ...)

        text(h$mids, h$counts, labels=h$counts, adj=c(0.5, -0.5), cex=.8)


    }else if(type=="mean"){

    ## histogram of mean coverage of subjects across loci
    cov_sub <- colMeans(coverage, na.rm=T)
    cov_sub_hist <- as.data.frame(cov_sub)

    h <- hist(cov_sub_hist[, 1],
              xlab="Mean coverage of subjects",
              main="Histogram of Mean Coverage of Subjects", ...)

    text(h$mids, h$counts, labels=h$counts, adj=c(0.5, -0.5), cex=.8)
    }

}
