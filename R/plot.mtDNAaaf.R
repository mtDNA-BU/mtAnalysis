#' plot.mtDNAaaf function
#'
#' Produce a scatter plot of a mtDNAaaf object
#'
#' @param x numeric matrix (16569 x N ). It contains subject ID as the column names,
#' and the AAF of all 16569 mtDNA loci for each subject.
#' @param col color to be used in the plot.
#' @param pch plotting ‘character’, i.e., symbol to use.
#' @param cex expansion factor for symbols used in the plot.
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param ... arguments to be passed to plot function.
#' @import graphics
#' @export
#' @examples
#'
#'
#' #aaf = mtAAF()
#' #plot (aaf)
#'
plot.mtDNAaaf <- function(x,  col = "blue", pch='.', cex=0.2, xlab="", ylab="", ...) {

  # faster approach - do not display zeros first:
  non.zeros <- which( x != 0 )
  plot(x = rep(c(1:16569),dim(x)[2])[non.zeros],
       y = as.vector(x[non.zeros]),
       col = col,
       pch = pch, cex = cex,
       xlab = xlab, ylab = ylab, ...)

  # add zeros
  points(x = 1:16569,
       y = rep(0, 16569),
       col = col,
       pch = pch, cex = cex,
       xlab = xlab, ylab = ylab, ...)

}
