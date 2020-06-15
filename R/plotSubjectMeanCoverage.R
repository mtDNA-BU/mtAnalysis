#' plotSubjectMeanCoverage function
#'
#' Produce a histogram of Mean Coverage Across Subjects
#'
#' @param aaf a numeric matrix (16569 x N). It contains subject ID as the row names, and the AAF of
#' all 16569 mtDNA loci for each subject.
#' @param coverage a numeric matrix (16569 x N) containing the subject ID as the row names, and the reads
#' coverage of the 16569 mitochondrial DNA loci for each subject.
#' @param coverage.qc a number(default is 100) of threshold for the coverage.
#' If the coverage<coverage.qc, the allele call at that locus of the subject will not be used.
#' @param thre.lower a number(default is 0.03) of lower bound of the threshold defining heteroplasmic
#' and homoplasmic variations
#' @param thre.upper a number(default is 0.97) of upper bound of the threshold defining heteroplasmic
#' and homoplasmic variations
#' @param loci one of: 1. a vector(default is c(1:16569)) of mitochondrial DNA loci to specify which loci
#' should be used to identify the variations and annotate, 2. a character string for the regions
#' (e.g. “coding”, “tRNA”, “RNR1”, "RNR2", …).
#' @import ggplot2
#' @export
#' @examples
#'
#'
#' #plotSubjectMeanCoverage (coverage)
#' #plotSubjectMeanCoverage (coverage, loci="tRNA")
#'
plotSubjectMeanCoverage <- function(aaf, coverage,
                                    coverage.qc=100,
                                    thre.lower=0.03,thre.upper=0.97,
                                    loci=c(1:16569)) {


  # give warning message and stop if the ncol of aaf, allele, freq and coverage are not 16569
  if(dim(aaf)[1]!=16569){
    stop("the coverage, allele, frequency and AAF should have 16569 loci (columns)")
  }

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


  # get the submatrices of aaf, allele, freq and coverage based on loci input
  aaf <- aaf[loci,]

  coverage <- coverage[loci,]
  rownames(coverage) <- as.character(loci)


  # categorize by AAF: AAF < thre.lower means no variation(coded as 0)
  # thre.lower <= AAF <= thre.upper means heteroplasmic variation(coded as 1)
  # thre.upper<AAF means homoplastic variation(coded as 2)

  aaf_cat <- matrix(0, nrow=nrow(aaf), ncol=ncol(aaf))
  aaf_cat[ aaf>=thre.lower & aaf<=thre.upper ] <- 1
  aaf_cat[ aaf> thre.upper ]    <- 2
  rownames(aaf_cat) <- rownames(aaf)
  colnames(aaf_cat) <- colnames(aaf)

  aaf_cat[ coverage < coverage.qc ]   <- NA

  # These are constants included in the package.
  # The reads of these loci are not reliable, so these loci should be removed
  loci_removed<-c(301, 302, 310, 316, 3107, 16182)
  loci_removed<-as.character(loci_removed)
  loci_left<-setdiff(rownames(aaf),loci_removed)
  aaf_cat <- aaf_cat[loci_left,]

  # Only include mutations loci: at least one subject has mutation at that locus
  mutation_collect <- aaf_cat
  # n_mutation_loci is the number of variations of each locus
  # find the loci which have at least 1 mutation (homo/heter) by n_mutation>0
  n_mutation<-rowSums(aaf_cat>0, na.rm=T)

  # only include the mutation(heter/homo) loci
  mutation_collect <- mutation_collect[ n_mutation > 0, ]

  # calculate the number of heteroplasmic variations for each subject
  heter_burden  <- colSums(mutation_collect == 1, na.rm=T)

  # calculate the number of heteroplasmic mutations for each locus
  heter_loci  <- rowSums(aaf_cat==1, na.rm=T)

  # calculate the number of homoplasmic variations for each subject
  homo_burden  <- colSums(mutation_collect==2, na.rm=T)

  # calculate the number of homoplasmic mutations for each locus
  homo_loci  <- rowSums(aaf_cat==2, na.rm=T)


  # generate histograms of heter and homo burden of across subjects and across mtDNA loci
  # histogram of heteroplasmic burden across subjects
  heter_burden_hist <- as.data.frame(heter_burden)
  p_heter_hist<-ggplot(heter_burden_hist, aes(x=heter_burden_hist[,1])) + geom_histogram(bins=100)
  p_heter_hist<-p_heter_hist+ theme_bw()
  p_heter_hist<-p_heter_hist+xlab("Histogram of heteroplasmic burden score")
  p_heter_hist<-p_heter_hist+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
  p_heter_hist<-p_heter_hist+coord_cartesian(xlim = c(0, 50))
  p_heter_hist <- p_heter_hist+stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5,bins=100)
  p_heter_hist

    # histogram of counts of heteroplasmic variations across mtDNA loci
    #heter_loci_hist<-as.data.frame(heter_loci)
    #p_heter_hist_loci<-ggplot(heter_loci_hist, aes(x=heter_loci_hist[,1])) + geom_histogram(bins=50)
    #p_heter_hist_loci<-p_heter_hist_loci+ theme_bw()
    #p_heter_hist_loci<-p_heter_hist_loci+xlab("Histogram of heteroplasmic variations of mtDNA loci")
    #p_heter_hist_loci<-p_heter_hist_loci+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
    #p_heter_hist_loci<-p_heter_hist_loci+coord_cartesian(xlim = c(0, 50))
    #heter_hist_loci<-p_heter_hist_loci+stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5,bins=50)
    #plot(p_heter_hist_loci)

    # histogram of homoplasmic burden across subjects
    #homo_burden_hist<-as.data.frame(homo_burden)
    #p_homo_hist<-ggplot(homo_burden_hist, aes(x=homo_burden_hist[,1])) + geom_histogram(bins=100)
    #p_homo_hist<-p_homo_hist+ theme_bw()
    #p_homo_hist<-p_homo_hist+xlab("Histogram of homoplasmic burden score")
    #p_homo_hist<-p_homo_hist+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
    #p_homo_hist<-p_homo_hist+coord_cartesian(xlim = c(0, 100))
    #homo_hist<-p_homo_hist+stat_bin(aes(y=..count.., label=..count..), geom="text", size=2.5,vjust=-1,bins=100)
    #plot(p_homo_hist)

    # histogram of counts of homoplasmic variations across mtDNA loci
    #homo_loci_hist<-as.data.frame(homo_loci)
    #p_homo_hist_loci<-ggplot(homo_loci_hist, aes(x=homo_loci_hist[,1])) + geom_histogram(bins=1000)
    #p_homo_hist_loci<-p_homo_hist_loci+ theme_bw()
    #p_homo_hist_loci<-p_homo_hist_loci+xlab("Histogram of homoplasmic variations of mtDNA loci")
    #p_homo_hist_loci<-p_homo_hist_loci+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
    #p_homo_hist_loci<-p_homo_hist_loci+coord_cartesian(xlim = c(0, 100))
    #homo_hist_loci<-p_homo_hist_loci+stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5,bins=1000)
    #plot(p_homo_hist_loci)



}
