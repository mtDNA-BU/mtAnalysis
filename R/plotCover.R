#' plotCover function
#'
#' Produce a Diagnostic Scatter Plot to Visualize Mean Coverage Across Subjects
#' for Each Locus
#'
#' @param x a numeric matrix (16569 x N).
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains the reads coverage of the 16569 mtDNA loci for each
#' subject. The matrix must contain the subject ID as the column names.
#' @param loci one of the following to specify mtDNA loci:
#' 1. a numeric vector (default is c(1:16569)) of mitochondrial DNA loci,
#' 2. a character string for the regions(e.g. "coding", "tRNA", "Dloop", …).
#' @param col color to be used in the plot.
#' @param pch plotting "character", i.e., symbol to use.
#' @param cex expansion factor for symbols used in the plot.
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param main a title for the plot.
#' @param type a string: "median" to calculate median coverage at individual
#' level, "mean" to calculate mean coverage at individual level
#' @param ... arguments to be passed to plot function.
#' @import graphics
#' @export
#' @examples
#'
#'\dontrun{
#' ## Read input data
#' coverage_file <- "coverage.csv"
#' coverage <- as.matrix(read.csv(file=coverage_file, sep=",", header=FALSE))
#'
#' plotCover(coverage, loci="coding")
#'
#'}
#'
plotCover <- function(x, loci=seq_len(.mtLength), col="blue", pch='.', cex=0.2,
                      xlab="mtloci",
                      ylab="", main="Mean Coverage for all mtDNA loci",
                      type="median",
                      ...) {

    if(length(type)!=1){
        stop("type must be a string of median or mean")
    }else if(type!="median" & type!="mean"){
        stop("type must be a string of median or mean")
    }

    if(dim(x)[1] != .mtLength)
        stop("the coverage should have 16569 loci (rows)")

    if(is.numeric(loci) & !all(loci %in% seq_len(.mtLength))){
        stop("loci should be a subset of 1:16569")
    }


    if (is.character(loci) ){
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

    if (is.character(x) )stop("coverage should be a numeric matrix")


    if(type=="median"){
        ## scatter plot of the median coverage across loci
        plot(x=loci,
             y=apply(x[loci, ], 2, FUN = median),
             col=col,
             pch=pch,
             cex=cex,
             xlab=xlab,
             ylab=ylab,
             main=main, ...)

    }else if(type=="mean"){
    ## scatter plot of the mean coverage across loci
    plot(x=loci,
         y=rowMeans(x[loci, ], na.rm=T),
         col=col,
         pch=pch,
         cex=cex,
         xlab=xlab,
         ylab=ylab,
         main=main, ...)

    }
}
