#' plotLociMeanCoverage function
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
#' #plotLociMeanCoverage (coverage)
#' #plotLociMeanCoverage (coverage, loci="tRNA")
#'
plotLociMeanCoverage <- function(x, loci=c(1:16569),
                             col = "blue", pch='.', cex=0.2, xlab="mtloci", ylab="",
                             main = "Mean Coverage across all mtDNA loci",
                             ...) {

  # give warning message and stop if the specified loci is not contained in 1:16569
  if(is.numeric(loci) & !all(loci %in% (1:16569))){
    stop("loci should be a subset of 1:16569")
  }


  #assign loci value when the specified loci are strings (e.g. coding)
  if (is.character(loci) ){
    if( loci=="coding"){
      loci<-sort(unique(c(3307:4262, 4470:5511,5904:7445,7586:8269,8366:8572,8527:9207
                          ,9207:9990 ,10059:10404, 10470:12137, 12337:14148, 14149:14673
                          , 14747:15887)))
    }else if(loci=="tRNA"){
      loci<-sort(unique(c(577:647, 1602:1670, 3230:3304, 4263:4331, 4329:4400, 4402:4469, 5512:5579
                          , 5587:5655, 5657:5729, 5761:5826, 5826:5891, 7446:7514, 7518:7585
                          , 8295:8364, 9991:10058, 10405:10469, 12138:12206, 12207:12265
                          , 12266:12336, 14674:14742, 15888:15953, 15956:16023)))
    }else if( loci=="RNR1"){
      loci<-sort(unique(c(648:1601)))
    }else if( loci=="RNR2"){
      loci<-sort(unique(c(1671:3229)))
    }
  }

  if (is.character(x) )stop("coverage should be a numeric matrix")


  # scatter plot of the mean coverage across loci
  plot(x = rowMeans(x[loci,], na.rm=T),
       y = loci,
       col = col,
       pch = pch,
       cex = cex,
       xlab = xlab,
       ylab = ylab,
       main = main, ...)


}
