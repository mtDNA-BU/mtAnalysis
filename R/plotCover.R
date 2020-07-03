#' plotCover function
#'
#' Produce a scatter plot of Mean Coverage Across Loci
#'
#' @param x a numeric matrix (16569 x N) containing the subject ID as the column names,
#' and the reads coverage of the 16569 mitochondrial DNA loci for each subject.
#' @param loci one of: 1. a vector(default is c(1:16569)) of mitochondrial DNA loci to specify which loci
#' should be used to identify the variations and annotate, 2. a character string for the regions
#' (e.g. “coding”, “tRNA”, “RNR1”, "RNR2", …).
#' @param col color to be used in the plot.
#' @param pch plotting ‘character’, i.e., symbol to use.
#' @param cex expansion factor for symbols used in the plot.
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param main a title for the plot.
#' @param ... arguments to be passed to plot function.
#' @import graphics
#' @export
#' @examples
#'
#'
#' #plotCover (coverage)
#' #plotCover (coverage, loci="tRNA")
#'
plotCover <- function(x, loci=c(1 : .mtLength),
                             col = "blue", pch='.', cex=0.2, xlab="mtloci", ylab="",
                             main = "Mean Coverage across all mtDNA loci",
                             ...) {

  # give warning message and stop if the specified loci is not contained in 1:16569
  if(is.numeric(loci) & !all(loci %in% (1 : .mtLength))){
    stop("loci should be a subset of 1:16569")
  }


  #assign loci value when the specified loci are strings (e.g. coding)
  if (is.character(loci) ){
    if( loci=="coding"){
      loci <- .loci.coding
    }else if(loci=="tRNA"){
      loci <- .loci.tRNA
    }else if( loci=="RNR1"){
      loci <- .loci.RNR1
    }else if( loci=="RNR2"){
      loci <- .loci.RNR2
    }else {
      stop("loci name must be one of: coding, tRNA, RNR1, RNR2")
    }
  }

  if (is.character(x) )stop("coverage should be a numeric matrix")


  # scatter plot of the mean coverage across loci
  plot(x = loci,
       y = rowMeans(x[loci,], na.rm=T),
       col = col,
       pch = pch,
       cex = cex,
       xlab = xlab,
       ylab = ylab,
       main = main, ...)


}
