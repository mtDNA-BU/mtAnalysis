#' histSampCov function
#'
#' Produce a Histogram of Mean Coverage Across Subjects
#'
#' @param coverage a numeric matrix (16569 x N) containing the subject ID as the row names, and the reads
#' coverage of the 16569 mitochondrial DNA loci for each subject.
#' @param loci one of: 1. a vector(default is c(1:16569)) of mitochondrial DNA loci to specify which loci
#' should be used to identify the variations and annotate, 2. a character string for the regions
#' (e.g. “coding”, “tRNA”, “RNR1”, "RNR2", …).
#' @import ggplot2
#' @export
#' @examples
#'
#'
#' #histSampCov (coverage)
#' #histSampCov (coverage, loci="tRNA")
#'
histSampCov <- function(coverage, loci=c(1:16569)) {


  # give warning message and stop if the ncol coverage are not 16569
  if(dim(coverage)[1]!=16569){
    stop("the coverage should have 16569 loci (columns)")
  }

  # give warning message and stop if the specified loci is not contained in 1:16569
  if(is.numeric(loci) & !all(loci %in% (1:16569))){
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


  # get the submatrices, allele, freq and coverage based on loci input
  coverage <- coverage[loci,]
  rownames(coverage) <- as.character(loci)


  # histogram of mean coverage of subjects across loci
  cov_sub<-colMeans(coverage , na.rm=T)
  cov_sub_hist <- as.data.frame(cov_sub)

  h <- hist(cov_sub_hist[,1], breaks = 50,
            xlab ="Mean coverage of subjects" ,
            main = "Histogram of Mean Coverage Across Subjects" )
  #rug( cov_sub_hist[,1] )
  text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5), cex = .8)

  # p_cov_hist <- ggplot(cov_sub_hist, aes(x=cov_sub_hist[,1])) +
  #   #geom_histogram(bins=100) +
  #   geom_bar(stat="bin") +
  #   geom_text(aes(label=len), vjust=-0.3, size=3.5)
  #   theme_bw() +
  #   xlab("Mean coverage of subjects") #+
  #   #stat_bin(aes(y=stat(count), label=stat(count)), geom="text", vjust=-.5,bins=50)
  # p_cov_hist

}
