#' mtSummary function
#' Identification of mtDNA variations, output summary statistics,
#' annotation of heteroplasmic and/or homoplasmic variations.
#' @param aaf a numeric matrix (N x 16569). It contains subject ID as the row names, and the AAF of
#' all 16569 mtDNA loci for each subject.
#' @param allele a data frame (N x 16569) provided by the user. This data frame contains
#' N subjects with mtDNA sequencing data of 16569 loci. The data frame contains subject ID as
#' the row names, and the allele calls of all mtDNA loci for each subject as the columns.
#' “/” is used to delimited different allele calls in a locus.
#' @param freq a data frame (N x 16569) provided by the user. This data frame contains the N subjects
#' with mtDNA sequencing data of 16569 loci. The data frame contains subject ID as the row names, and
#' the allele fractions of the called alleles for each subject as the columns. “/” is used to delimited
#' the allele fractions.
#' @param coverage a numeric matrix (N x 16569) containing the subject ID as the row names, and the reads
#' coverage of the 16569 mitochondrial DNA loci for each subject.
#' @param coverage.qc a number(default is 100) of threshold for the coverage.
#' If the coverage<coverage.qc, the allele call at that locus of the subject will not be used.
#' @param thre.lower a number(default is 0.03) of lower bound of the threshold defining heteroplasmic
#' and homoplasmic variations
#' @param thre.upper a number(default is 0.97) of upper bound of the threshold defining heteroplasmic
#' and homoplasmic variations
#' @param loci one of: 1. a vector(default is c(1:16569)) of mitochondrial DNA loci to specify which loci
#' should be used to identify the variations and annotate, 2. a character string for the regions
#' (e.g. “coding”, “tRNA”, “Dloop”, …).
#' @param type a number of indicator choosing to output annotation to all variations, heteroplasmic
#' variations, or homoplasmic variations. 0 returns annotation to all variations (default), 1 returns
#' annotation to heteroplasmic variations and 2 returns annotation to homoplasmic variations.
#' @param coverSummary logical(default is True). A user can specify to output scatter plot of mean
#' coverage across mtDNA loci and summary of coverage at each mtDNA loci (across all participants) and
#' for each individual (across all mtDNA loci).
#' @param varHist logical(default is True). A user can specify to output histograms to visualize the
#' heteroplasmic and homoplasmic burden across participants and mtDNA loci.
#' @param annot.select types of annotation scores based on user's choice
#' @param path the path of the output annotation file. If not provided, the annotation file will
#' output to the current working directory
#' @param study A string of study names. Default is Study.
#'
#' @return 1. Summary frequency: descriptive statistics of sequencing coverage across individual
#' and across mtDNA loci; the total number of mtDNA loci with variations and the number of
#' heteroplasmic/homoplasmic mtDNA loci in the study sample; the min/Q1/mean/median/Q3/largest number
#' of homoplasmy/heteroplasmy carried by an individual. 2. Plots: a scatter plot of the median coverage
#' across mtDNA loci; histograms to visualize the heteroplasmic and homoplasmic burden across
#' participants and mtDNA loci based on user’s choice. 3.	Summary annotation for all of the
#' heteroplasmic and/or homoplasmic variations observed in the study data based on user’s choice.
#' @import graphics ggplot2 utils
#' @export
#' @examples
#'
#' #mtSummary(allele, freq, coverage, thre.lower, thre.upper, loci, type, coverSummary, varHist)
mtSummary<-function(aaf, allele,freq,coverage,coverage.qc=100,thre.lower=0.03,thre.upper=0.97,loci=c(1:16569),type=0,coverSummary=T,varHist=T,
                    annot.select=c("Pos","ref","Gene","TypeMutation","MissensMutation","CodonPosition","ProteinDomain","dbSNP_150_id","PolyPhen2","PolyPhen2_score","SIFT","SIFT_score", "CADD","CADD_score","CADD_phred_score"),
                    path="",study="Study"){
  # give warning message and stop if the aaf, allele, freq and coverage do not have the same dimension
  if(!all(sapply(list(dim(aaf),dim(allele),dim(freq),dim(coverage)), function(x) x == dim(coverage)))){
    stop("the coverage, allele and frequency should have same dimension")
  }

  # give warning message and stop if the length of loci is larger than 16569
  # the total number of loci of mitochondrial DNA
  if(length(loci)>16569){
    stop("loci should have no more than 16569")
  }

  # transpose the aaf, allele, freq and coverage dataset
  aaf<-t(aaf)
  allele<-t(allele)
  freq<-t(freq)
  coverage<-t(coverage)

  # get the submatrices of aaf, allele, freq and coverage based on loci input
  aaf<-aaf[loci,]
  allele<-allele[loci,]
  freq<-freq[loci,]
  coverage<-coverage[loci,]

  # create a object of output of summary
  mt_summary_obj<-NULL

  if(coverSummary){
    coverage_mean_loci<-apply(coverage, 1, mean)
    # output the summary of mean coverage of each locus
    mt_summary_obj$coverLoci<-summary(coverage_mean_loci)
    #length(coverage_mean_loci)
    coverage_mean_subjects<-apply(coverage, 2, mean)
    # output the summary of mean coverage of each subject
    mt_summary_obj$coverSubjects<-summary(coverage_mean_subjects)
    #length(coverage_mean_subjects)
    # scatter plot of the mean coverage across loci
    y_coverage_scatter<-rowMeans(coverage,na.rm=T)
    x_coverage_scatter<-seq(1:dim(coverage)[1])
    plot(x=x_coverage_scatter, y=y_coverage_scatter,col = "blue",pch = ".",cex=0.2,xlab="mtloci",ylab="",main="mean coverage across all mtDNA loci")
    # scatter plot of the mean coverage across subjects
    # y_coverage_scatter<-colMeans(coverage,na.rm=T)
    # x_coverage_scatter<-seq(1:dim(coverage)[2])
    # plot(x=x_coverage_scatter, y=y_coverage_scatter,col = "blue",pch = ".",cex=0.2,xlab="mtloci",ylab="",main="mean coverage across all subjects")
  }

  # categorize by AAF: AAF<thre.lower means no variation(coded as 0)
  # thre.lower<=AAF<=thre.upper means heteroplasmic variation(coded as 1)
  # thre.upper<AAF means homoplastic variation(coded as 2)
  aaf_cat<-aaf
  aaf_cat[ aaf< thre.lower ]    <- 0
  aaf_cat[ aaf>=thre.lower & aaf<=thre.upper ] <- 1
  aaf_cat[ aaf> thre.upper ]    <- 2

  # if coverage<coverage.qc, means unreliable read, so assign aaf_cat to NA
  aaf_cat[ coverage < coverage.qc ]       <- NA

  # These are constants included in the package.
  # The reads of these loci are not reliable, so these loci should be removed
  loci_removed<-c(301,302,310,316,3107,16182)
  aaf_cat      <- aaf_cat[-loci_removed,]

  # Only include mutations loci: at least one subject has mutation at that locus
  mutation_collect<-aaf_cat
  # n_mutation_loci is the number of variations of each locus
  # find the loci which have at least 1 mutation (homo/heter) by n_mutation>0
  n_mutation<-rowSums(aaf_cat>0, na.rm=T)
  # length(n_mutation)
  # only include the mutation(heter/homo) loci
  mutation_collect<-mutation_collect[n_mutation>0,]
  #dim(mutation_collect)
  # loci of variations(heter/homo)
  loci_var<-as.numeric(rownames(mutation_collect))
  mt_summary_obj$loci_var<-loci_var
  n_mutation<-n_mutation[n_mutation>0]
  # length(n_mutation)
  # [1] 4193


  # output summary statistics of hereoplasmic variations
  # calculate the number of heteroplasmic variations for each subject
  heter_burden  <- colSums(mutation_collect==1, na.rm=T)
  # output the summary of heteroplasmic burden
  mt_summary_obj$heter_burden_sum<-summary(heter_burden)
  # calculate the number of heteroplasmic mutations for each locus
  heter_loci  <- rowSums(aaf_cat==1, na.rm=T)
  # output the summary of heteroplasmic variations across loci
  mt_summary_obj$heter_loci_sum<-summary(heter_loci)
  # loci with heteroplasmic variations
  loci_heter<-as.numeric(names(heter_loci))[heter_loci>0]
  mt_summary_obj$loci_heter<-loci_heter
  # calculate the total number of heteroplasmic variations
  mt_summary_obj$heter_total<-sum(heter_burden)

  # output summary statistics of homoplasmic variations
  # calculate the number of homoplasmic variations for each subject
  homo_burden  <- colSums(mutation_collect==2, na.rm=T)
  # output the summary homoplasmic burden
  mt_summary_obj$homo_burden<-summary(homo_burden)
  # calculate the number of homoplasmic mutations for each locus
  homo_loci  <- rowSums(aaf_cat==2, na.rm=T)
  # output the summay of homoplasmic variations across loci
  mt_summary_obj$homo_loci_sum<-summary(homo_loci)
  # loci with homoplasmic variations
  loci_homo<-as.numeric(names(homo_loci))[homo_loci>0]
  mt_summary_obj$loci_homo<-loci_homo
  # calculate the total number of homoplasmic variations
  mt_summary_obj$homo_total<-sum(homo_burden)

  # generate histograms of heter/homo burden of acorss subjects and across mtDNA loci
  if(varHist){
    # histogram of heteroplasmic burden across subjects
    heter_burden_hist<-as.data.frame(heter_burden)
    p_heter_hist<-ggplot(heter_burden_hist, aes(x=heter_burden_hist[,1])) + geom_histogram(bins=100)
    p_heter_hist<-p_heter_hist+ theme_bw()
    p_heter_hist<-p_heter_hist+xlab("Histogram of heteroplasmic burden score")
    p_heter_hist<-p_heter_hist+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
    p_heter_hist<-p_heter_hist+coord_cartesian(xlim = c(0, 50))
    p_heter_hist<-p_heter_hist+stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5,bins=100)
    plot(p_heter_hist)

    # histogram of counts of heteroplasmic variations across mtDNA loci
    heter_loci_hist<-as.data.frame(heter_loci)
    p_heter_hist_loci<-ggplot(heter_loci_hist, aes(x=heter_loci_hist[,1])) + geom_histogram(bins=1000)
    p_heter_hist_loci<-p_heter_hist_loci+ theme_bw()
    p_heter_hist_loci<-p_heter_hist_loci+xlab("Histogram of heteroplasmic variations across mtDNA loci")
    p_heter_hist_loci<-p_heter_hist_loci+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
    p_heter_hist_loci<-p_heter_hist_loci+coord_cartesian(xlim = c(0, 50))
    p_heter_hist_loci<-p_heter_hist_loci+stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5,bins=1000)
    plot(p_heter_hist_loci)

    # histogram of homoplasmic burden across subjects
    homo_burden_hist<-as.data.frame(homo_burden)
    p_homo_hist<-ggplot(homo_burden_hist, aes(x=homo_burden_hist[,1])) + geom_histogram(bins=100)
    p_homo_hist<-p_homo_hist+ theme_bw()
    p_homo_hist<-p_homo_hist+xlab("Histogram of homoplasmic burden score")
    p_homo_hist<-p_homo_hist+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
    p_homo_hist<-p_homo_hist+coord_cartesian(xlim = c(0, 100))
    p_homo_hist<-p_homo_hist+stat_bin(aes(y=..count.., label=..count..), geom="text", size=2.5,vjust=-1,bins=100)
    plot(p_homo_hist)

    # histogram of counts of homoplasmic variations across mtDNA loci
    homo_loci_hist<-as.data.frame(homo_loci)
    p_homo_hist_loci<-ggplot(homo_loci_hist, aes(x=homo_loci_hist[,1])) + geom_histogram(bins=1000)
    p_homo_hist_loci<-p_homo_hist_loci+ theme_bw()
    p_homo_hist_loci<-p_homo_hist_loci+xlab("Histogram of homoplasmic variations across mtDNA loci")
    p_homo_hist_loci<-p_homo_hist_loci+theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=15),axis.title=element_text(size=15,face="bold"))
    p_homo_hist_loci<-p_homo_hist_loci+coord_cartesian(xlim = c(0, 100))
    p_homo_hist_loci<-p_homo_hist_loci+stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5,bins=1000)
    plot(p_homo_hist_loci)
  }

  # annotation for all heter and/or homo variations at each mutation loci based on user's choice
  # if type==0, annotate both heter/homo variations
  if(type==0){
    # allele2 is the submatrix of allele which only includes variation loci
    allele2<-allele[loci_var,]
    # freq2 is the submatrix of freq which only includes variation loci
    freq2<-freq[loci_var,]

    allele_both<-rep(NA,dim(mutation_collect)[1])
    # identify all the alleles(including ref allele and alternative alleles)
    # which the frequency falls in interval [thre.lower, 1] for each variaation locus
    for(i in 1:dim(mutation_collect)[1]){
      var_index<-(mutation_collect[i,]>0)
      # identify the alleles at the variation locus
      allele_all<-allele2[i,]
      allele_var<-allele_all[var_index]
      allele_var <- allele_var[!is.na(allele_var)]

      # identify the corresponding freqencies of the alleles at the variation loci
      freq_all<-freq2[i,]
      freq_var<-freq_all[var_index]
      freq_var <- freq_var[!is.na(freq_var)]

      allele_all_var<-unlist(strsplit(allele_var,split="/"))
      freq_all_var<-as.numeric(unlist(strsplit(freq_var ,split="/")))
      # only includes alleles which have frequency within the [thre.lower, 1] interval
      allele_all_var<-allele_all_var[freq_all_var>=thre.lower]
      # combine with the reference allele
      allele_all_var<-c(.mtRef[as.numeric(rownames(mutation_collect)[i])],allele_all_var)

      allele_all_var<-allele_all_var[!duplicated(allele_all_var)]
      allele_all_var<-paste(allele_all_var,collapse = '/')
      allele_both[i]<-allele_all_var
    }

    # choose the type of scores to be annotated
    # annot.select: they are provided by users with default
    # annot.select <- c("Pos","ref","Gene","TypeMutation","MissensMutation","CodonPosition","ProteinDomain","dbSNP_150_id","PolyPhen2","PolyPhen2_score","SIFT","SIFT_score", "CADD","CADD_score","CADD_phred_score")
    head         <- c("mtID","ref","allele_var","n_var","n_heter","n_homo","mut_allele",annot.select)
    # output the annotation as a .csv file to the path users provided
    write (head, file= paste0(path,study,"_annotation_variation.csv"),sep=",",ncolumns=length(head), append=T)

    # this loop used to annote all the mutations at each heter loci with the scores chosen(annot.select)
    for (i in 1:dim(mutation_collect)[1]) {
      # identify heter locus
      pos       <- loci_var[i]

      # assign the ref allele at that position
      ref_allele <- .mtRef[pos]

      # assign all the heter alleles at that locus
      x    <- allele_both[i]
      x2   <- unlist(strsplit(x,split="/"))

      # if statement to check if length(x2)>1
      if(length(x2)>1){

      # find the positions that the alleles which are not ref allele
      pos.not  <- which(!(x2==ref_allele))
      # this loop is used to annotate all the alternative alleles(with freq>=thre.lower) at that locus
      # according to the AA, CC, GG, TT annotation files
      for (k in 1:length(pos.not)){
        point    <- x2[pos.not[k]]
        if (which (point == c("A","C","G","T"))==1) {
          score   <- AA[AA$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,n_mutation[i],heter_loci[i],homo_loci[i],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_variation.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
        else if (which (point == c("A","C","G","T"))==2) {
          score <- CC[CC$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,n_mutation[i],heter_loci[i],homo_loci[i],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_variation.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)

        }
        else if (which (point == c("A","C","G","T"))==3) {
          score <- GG[GG$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,n_mutation[i],heter_loci[i],homo_loci[i],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_variation.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
        else if (which (point == c("A","C","G","T"))==4) {
          score   <- TT[TT$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,n_mutation[i],heter_loci[i],homo_loci[i],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_variation.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
      }
      }
    }
  }

  # if type==1, annotate heter variations
  else if(type==1){
    # allele2 is a submatrix of allele which only includes heteroplasmic loci
    allele2<-allele[loci_heter,]
    # freq2 is a submatrix of freq which only includes heteroplasmic loci
    freq2<-freq[loci_heter,]
    # heter_collect is a submatrix of mutation_collect which only includes heteroplasmic loci
    heter_collect<-mutation_collect[as.numeric(rownames(mutation_collect))%in%loci_heter,]

    allele_heter<-rep(NA,length(loci_heter))
    # identify all the alleles(including ref allele and alternative alleles)
    # which the frequency falls in interval [thre.lower, thre.upper] for each variation locus
    for(i in 1:dim(heter_collect)[1]){
      heter_index<-(heter_collect[i,]==1)
      # identify the alleles at the variation locus
      allele_all<-allele2[i,]
      allele_var<-allele_all[heter_index]
      allele_var <- allele_var[!is.na(allele_var)]

      # identify the corresponding freqencies of the alleles at the variation loci
      freq_all<-freq2[i,]
      freq_var<-freq_all[heter_index]
      freq_var <- freq_var[!is.na(freq_var)]

      allele_all_var<-unlist(strsplit(allele_var,split="/"))
      freq_all_var<-as.numeric(unlist(strsplit(freq_var ,split="/")))
      # only includes alleles which have frequency within the [thre.lower, thre.upper] interval
      allele_all_var<-allele_all_var[freq_all_var>=thre.lower & freq_all_var<=thre.upper]
      # combine with the reference allele
      allele_all_var<-c(.mtRef[as.numeric(rownames(heter_collect)[i])],allele_all_var)

      allele_all_var<-allele_all_var[!duplicated(allele_all_var)]
      allele_all_var<-paste(allele_all_var,collapse = '/')
      allele_heter[i]<-allele_all_var
    }

    # choose the type of scores to be annotated
    # annot.select: they are provided by users with default
    # annot.select <- c("Pos","ref","Gene","TypeMutation","MissensMutation","CodonPosition","ProteinDomain","dbSNP_150_id","PolyPhen2","PolyPhen2_score","SIFT","SIFT_score", "CADD","CADD_score","CADD_phred_score")
    head         <- c("mtID","ref","allele_heter","n_heter","mut_allele",annot.select)
    # output the annotation as a .csv file to the path users provided
    write (head, file= paste0(path,study,"_annotation_heteroplasmy.csv"),sep=",",ncolumns=length(head), append=T)

    # this loop used to annote all the heter mutations at each heter loci with the scores chosen(annot.select)
    for (i in 1:dim(heter_collect)[1]) {
      # identify heter locus
      pos       <- loci_heter[i]

      # assign the ref allele at that position
      ref_allele <- .mtRef[pos]

      # assign all the heter alleles at that locus
      x    <- allele_heter[i]
      x2   <- unlist(strsplit(x,split="/"))

      # if statement to check if length(x2)>1
      if(length(x2)>1){

      # find the positions that the alleles which are not ref allele
      pos.not  <- which(!(x2==ref_allele))
      # this loop is used to annotate all the heter alleles at that locus(excluding the  ref allele)
      # according to the AA, CC, GG, TT annotation files
      for (k in 1:length(pos.not)){
        point    <- x2[pos.not[k]]
        if (which (point == c("A","C","G","T"))==1) {
          score   <- AA[AA$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,heter_loci[loci_var==loci_heter[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_heteroplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
        else if (which (point == c("A","C","G","T"))==2) {
          score <- CC[CC$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,heter_loci[loci_var==loci_heter[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_heteroplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)

        }
        else if (which (point == c("A","C","G","T"))==3) {
          score <- GG[GG$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,heter_loci[loci_var==loci_heter[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_heteroplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
        else if (which (point == c("A","C","G","T"))==4) {
          score   <- TT[TT$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,heter_loci[loci_var==loci_heter[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_heteroplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
      }
      }
    }
  }

  # if type==2, annotate homo variations
  else if(type==2){
    # allele2 is a submatrix of allele which only includes homoplasmic loci
    allele2<-allele[loci_homo,]
    # freq2 is a submatrix of freq which only includes homoplasmic loci
    freq2<-freq[loci_homo,]
    # homo_collect is a submatrix of mutation_collect which only includes homoplasmic loci
    homo_collect<-mutation_collect[as.numeric(rownames(mutation_collect))%in%loci_homo,]

    allele_homo<-rep(NA,length(loci_homo))
    # identify all the alleles(including ref allele and alternative alleles)
    # which the frequency falls in interval (thre.upper, 1] for each variation locus
    for(i in 1:dim(homo_collect)[1]){
      homo_index<-(homo_collect[i,]==1)
      # identify the alleles at the variation locus
      allele_all<-allele2[i,]
      allele_var<-allele_all[homo_index]
      allele_var <- allele_var[!is.na(allele_var)]

      # identify the corresponding freqencies of the alleles at the variation loci
      freq_all<-freq2[i,]
      freq_var<-freq_all[homo_index]
      freq_var <- freq_var[!is.na(freq_var)]

      allele_all_var<-unlist(strsplit(allele_var,split="/"))
      freq_all_var<-as.numeric(unlist(strsplit(freq_var ,split="/")))
      # only includes alleles which have frequency within the [thre.lower, thre.upper] interval
      allele_all_var<-allele_all_var[freq_all_var>thre.upper]
      # combine with the reference allele
      allele_all_var<-c(.mtRef[as.numeric(rownames(homo_collect)[i])],allele_all_var)

      allele_all_var<-allele_all_var[!duplicated(allele_all_var)]
      allele_all_var<-paste(allele_all_var,collapse = '/')
      allele_homo[i]<-allele_all_var
    }

    # choose the type of scores to be annotated
    # annot.select: they are provided by users with default
    # annot.select <- c("Pos","ref","Gene","TypeMutation","MissensMutation","CodonPosition","ProteinDomain","dbSNP_150_id","PolyPhen2","PolyPhen2_score","SIFT","SIFT_score", "CADD","CADD_score","CADD_phred_score")
    head         <- c("mtID","ref","allele_homo","n_homo","mut_allele",annot.select)
    # output the annotation as a .csv file to the path users provided
    write (head, file= paste0(path,study,"_annotation_homoplasmy.csv"),sep=",",ncolumns=length(head), append=T)

    # this loop used to annote all the heter mutations at each heter loci with the scores chosen(annot.select)
    for (i in 1:dim(homo_collect)[1]) {
      # identify heter locus
      pos       <- loci_homo[i]

      # assign the ref allele at that position
      ref_allele <- .mtRef[pos]

      # assign all the heter alleles at that locus
      x    <- allele_homo[i]
      x2   <- unlist(strsplit(x,split="/"))

      # if statement to check if length(x2)>1
      if(length(x2)>1){

      # find the positions that the alleles which are not ref allele
      pos.not  <- which(!(x2==ref_allele))
      # this loop is used to annotate all the heter alleles at that locus(excluding the  ref allele)
      # according to the AA, CC, GG, TT annotation files
      for (k in 1:length(pos.not)){
        point    <- x2[pos.not[k]]
        if (which (point == c("A","C","G","T"))==1) {
          score   <- AA[AA$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,homo_loci[loci_var==loci_homo[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_homoplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
        else if (which (point == c("A","C","G","T"))==2) {
          score <- CC[CC$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,homo_loci[loci_var==loci_homo[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_homoplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)

        }
        else if (which (point == c("A","C","G","T"))==3) {
          score <- GG[GG$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,homo_loci[loci_var==loci_homo[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_homoplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
        else if (which (point == c("A","C","G","T"))==4) {
          score   <- TT[TT$Pos==pos,annot.select]
          all   <- cbind(pos,ref_allele,x,homo_loci[loci_var==loci_homo[i]],point ,score)
          write.table(all, file= paste0(path,study,"_annotation_homoplasmy.csv"),sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
      }
      }
    }
  }
  return(mt_summary_obj)
}










