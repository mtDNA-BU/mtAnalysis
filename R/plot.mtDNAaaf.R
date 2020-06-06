#' plot.mtDNAaaf function
#'
#' Plot mtDNAaaf object
#'
#' @param x numeric matrix (N x 16569). It contains subject ID as the row names,
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
#' #plot (aaf)
#'
plot.mtDNAaaf <- function(x,  col = "blue", pch='.', cex=0.2, xlab="", ylab="", ...) {

  #AAF_scatter_x <- rep(c(1:16569),dim(x)[2])
  #AAF_scatter_y<- as.vector( x )
  #plot(x=AAF_scatter_x, y=AAF_scatter_y,
  #     col = col,
  #     pch = pch, cex = cex,
  #     xlab = xlab, ylab = ylab, ...)

  # faster approach - do not display zeros:
  non.zeros <- which( x!=0 )
  plot(x = rep(c(1:16569),dim(x)[2])[non.zeros],
       y = as.vector(x[non.zeros]),
       col = col,
       pch = pch, cex = cex,
       xlab = xlab, ylab = ylab, ...)
  # add zeros
  plot(x = 1:16569, y = 0,
       add = TRUE,
       col = col,
       pch = pch, cex = cex,
       xlab = xlab, ylab = ylab, ...)

}
