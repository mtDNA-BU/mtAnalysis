#' mtSummary function
#' Identification of mtDNA variations, output summary statistics,
#' annotation of heteroplasmic and/or homoplasmic variations.
#' @param aaf a numeric matrix (16569 x N). Rows correspond to loci and columns correspond to subjects.
#' It contains subject ID as the column names, and the AAFs of all 16569 mtDNA loci for each subject.
#' It is generated from mtAAF function.
#' @param allele a data frame (16569 x N) provided by the user. Rows correspond to loci and columns
#' correspond to subjects. This data frame contains N subjects with mtDNA sequencing data of 16569 loci.
#' The data frame contains subject ID as the column names. “/” is used to delimited different allele calls
#' in a locus.
#' @param freq a data frame (16569 x N) provided by the user. Rows correspond to loci and columns correspond
#' to subjects. This data frame contains the N subjects with mtDNA sequencing data of 16569 loci. The data
#' frame contains subject ID as the column names. “/” is used to delimited the allele fractions.
#' @param coverage a numeric matrix (16569 x N). Rows correspond to loci and columns correspond to subjects.
#' It contains the subject ID as the column names, and the reads coverage of the 16569 mtDNA loci for each
#' subject.
#' @param coverage.qc a number (default is 100) of threshold for the coverage. If the coverage<coverage.qc,
#' the allele call at that locus of the participant will not be used.
#' @param thre.lower a number(default is 0.03) of lower bound of the threshold defining heteroplasmic and
#' homoplasmic variations.
#' @param thre.upper a number(default is 0.97) of upper bound of the threshold defining heteroplasmic and
#' homoplasmic variations.
#' @param loci one of the following to specify which loci should be used to identify the variations and
#' annotate: 1. a numeric vector (default is c(1:16569)) of mitochondrial DNA loci, 2. a character string
#' for the regions(e.g. “coding”, “tRNA”, “Dloop”, …).
#' @param type a character string of indicator choosing to output annotation to all variations,
#' heteroplasmic mutations, or homoplasmic variations. “both” returns annotation to all variations (default),
#'  “heter” returns annotation to heteroplasmic variations and “homo” returns annotation to homoplasmic
#'  variations.
#' @param coverSummary : logical (default is True). A user can specify to output summary of mean coverage
#' across mtDNA loci, and summary of mean coverage across all subjects.
#' @param varHist logical (default is True). A user can specify to output histograms to visualize the
#' heteroplasmic and homoplasmic burden of participants and number of mutations observed at each mtDNA locus.
#' @param annot.select types of annotation scores to be generated.
#' @param path the path to the directory of the output of annotation file. If not provided, the annotation
#' file will output to the current working directory.
#' @param study A string of study names. Default is “Study”.
#'
#' @return 1. A list containing summary of mean coverage across mtDNA loci, summary of mean coverage across
#' all subjects, loci of heteroplasmic and homoplasmic mutations, summary of heteroplasmic and homoplasmic
#' burden of subjects, summary of number of heteroplasmic and homoplasmic mutations of loci, and total
#' number of heteroplasmic and homoplasmic mutations. 2. Histograms of heteroplasmic and homoplasmic burden
#' of participants and number of mutations observed at mtDNA loci. 3. A .csv file containing annotated
#' alternative of heteroplasmic and/or homoplasmic variations observed in the data with corresponding mtDNA
#' loci positions.
#' @import graphics grDevices utils
#' @export
#' @examples
#'
#' #mtSummary(aaf, allele, freq, coverage, thre.lower, thre.upper, loci, type, coverSummary, varHist)
mtSummary<-function(aaf, allele, freq, coverage,
                    coverage.qc=100, thre.lower=0.03,thre.upper=0.97,
                    loci=c(1 : .mtLength),
                    type="both",
                    coverSummary=T, varHist=T,
                    annot.select=c("Pos","ref","Gene","TypeMutation","MissensMutation",
                                   "CodonPosition","ProteinDomain","dbSNP_150_id","PolyPhen2",
                                   "PolyPhen2_score","SIFT","SIFT_score", "CADD","CADD_score",
                                   "CADD_phred_score"),
                    path="./", study="Study"){

  if(!all(is.character(allele) , is.character(freq) , is.numeric(aaf) , is.numeric(coverage))){
    stop("the allele and freq shoud be character, aaf and coverage should be numeric")
  }

  if(!all(sapply(list(dim(aaf), dim(allele), dim(freq)), function(x) x == dim(coverage)))){
    stop("the coverage, allele, frequency and AAF should have same dimension")
  }

  if(dim(aaf)[1]!=.mtLength){
    stop("the coverage, allele, frequency and aaf should have 16569 loci (columns)")
  }

  if(!all(setequal(colnames(allele),colnames(freq)) , setequal(colnames(allele),colnames(aaf)), setequal(colnames(allele),colnames(coverage)))){
    stop("the coverage, allele, frequency and aaf should have same subject IDs (columns)")
  }


  if(!all(mode(loci)%in%c("numeric","character") , is.vector(loci))){
    stop("loci must be numeric vector or character")
  }

  if(is.character(loci)){
    if(!(all(loci%in%c("coding", "tRNA", "RNR1", "RNR2")) & length(loci)==1)){
    stop("loci must be one of coding, tRNA, RNR1, RNR2 if it is character")
    }
  }

  if(is.numeric(loci)){
    if(!all(loci %in% (1:.mtLength))){
    stop("loci should be a subset of 1:16569")
    }
  }

  if( ! dir.exists(path) ) stop( paste(" Output path", path,"does not exist.") )

  if(!all(is.numeric(coverage.qc) , length(coverage.qc)==1)){
    stop("coverage.qc must be numeric of length 1")
  }

  if(!all(is.numeric(c(thre.lower,thre.upper)) , length(thre.lower)==1 , length(thre.upper)==1)){
    stop("thre.lower and thre.upper must be numeric of length 1")
  }

  if(thre.lower>=thre.upper){
    stop("thre.lower must be smaller than thre.upper")
  }

  if(!all(length(type)==1 , type %in% c("both","heter","homo"))){
    stop("type must be one of both, heter and homo")
  }

  # Check if annot.select belongs to the colnames of annotation files
  if(!all(annot.select %in% names(.AA) , is.vector(annot.select))){
    stop("annot.select must be a subset of variables of annotation files")
  }

  subjectID <- colnames(allele)
  subjectID<-sort(subjectID)
  allele <- allele[ , subjectID ]
  freq <- freq[ , subjectID]
  coverage <-coverage[ , subjectID]
  aaf <-aaf[ , subjectID]

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

  aaf <- aaf[loci,]
  allele <- allele[loci,]
  rownames(allele) <- as.character(loci)

  freq <- freq[loci,]
  rownames(freq) <- as.character(loci)

  coverage <- coverage[loci,]
  rownames(coverage) <- as.character(loci)

  # create an object of output of summary
  mt_summary_obj <- NULL


  if(coverSummary){
    mt_summary_obj$coverLoci <-     summary( rowMeans(coverage, na.rm=T) )
    mt_summary_obj$coverSubjects <- summary( colMeans(coverage, na.rm=T) )
  }

  # categorize by AAF: AAF < thre.lower means no variation(coded as 0)
  # thre.lower <= AAF <= thre.upper means heteroplasmic variation(coded as 1)
  # thre.upper<AAF means homoplastic variation(coded as 2)

  aaf_cat <- matrix(0, nrow=nrow(aaf), ncol=ncol(aaf))
  aaf_cat[ aaf>=thre.lower & aaf<=thre.upper ] <- 1
  aaf_cat[ aaf> thre.upper ]    <- 2
  rownames(aaf_cat) <- rownames(aaf)
  colnames(aaf_cat) <- colnames(aaf)


  aaf_cat[ coverage < coverage.qc ]   <- NA

  # The reads of these loci are not reliable, so these loci should be removed
  loci_removed <- .mtLociNUMT   # c(301,302,310,316,3107,16182)
  loci_removed <- as.character(loci_removed)
  loci_left <- setdiff(rownames(aaf),loci_removed)
  aaf_cat      <- aaf_cat[loci_left,]

  # only include mutations loci: at least one subject has mutation at that locus
  mutation_collect <- aaf_cat

  n_mutation <- rowSums(aaf_cat>0, na.rm=T)

  # only include the mutation(heter/homo) loci
  mutation_collect<-mutation_collect[n_mutation>0,]

  # loci of variations(heter/homo)
  loci_var<-as.numeric(rownames(mutation_collect))
  mt_summary_obj$loci_var<-loci_var
  n_mutation <- n_mutation[n_mutation>0]

  # output summary statistics of hereoplasmic variations
  # calculate the number of heteroplasmic mutations for each subject
  heter_burden  <- colSums(aaf_cat==1, na.rm=T)
  mt_summary_obj$heter_burden_sum<-summary(heter_burden)

  # calculate the number of heteroplasmic mutations for each locus
  heter_loci  <- rowSums(aaf_cat==1, na.rm=T)
  mt_summary_obj$heter_loci_sum<-summary(heter_loci)

  # loci with heteroplasmic variations
  loci_heter<-as.numeric(names(heter_loci))[heter_loci>0]
  mt_summary_obj$loci_heter<-loci_heter

  # calculate the total number of heteroplasmic variations
  mt_summary_obj$heter_total<-sum(heter_burden)

  # output summary statistics of homoplasmic variations
  # calculate the number of homoplasmic variations for each subject
  homo_burden  <- colSums(aaf_cat==2, na.rm=T)
  mt_summary_obj$homo_burden_sum<-summary(homo_burden)

  # calculate the number of homoplasmic mutations for each locus
  homo_loci  <- rowSums(aaf_cat==2, na.rm=T)
  mt_summary_obj$homo_loci_sum<-summary(homo_loci)

  # loci with homoplasmic variations
  loci_homo<-as.numeric(names(homo_loci))[homo_loci>0]
  mt_summary_obj$loci_homo<-loci_homo

  # calculate the total number of homoplasmic variations
  mt_summary_obj$homo_total<-sum(homo_burden)

  # generate histograms of heter and homo burden of across subjects and across mtDNA loci
  if(varHist){

    #Open an output pdf file
    pdf(file = paste0(path,"/", study,"_mtHistograms.pdf"))

    # histogram of heteroplasmic burden across subjects
    h <- hist(heter_burden, breaks = 100,
              xlab ="Heteroplasmic burden score" ,
              main = "Histogram of heteroplasmic burden score" )

    h <- hist(heter_loci, breaks = 50,
              xlab ="heteroplasmic variations of mtDNA loci" ,
              main = "Histogram of heteroplasmic variations of mtDNA loci" )

    h <- hist(homo_burden, breaks = 100,
              xlab ="Homoplasmic burden score" ,
              main = "Histogram of homoplasmic burden score" )

    h <- hist(homo_loci, breaks = 600, xlim = c(0, 150),
              xlab ="Homoplasmic variations of mtDNA loci" ,
              main = "Histogram of homoplasmic variations of mtDNA loci" )

    dev.off()
  }

  # annotation for all heter and/or homo variations at each mutation loci based on user's choice
  # if type=="both", annotate both heter/homo variations
  if(type=="both"){

    allele2<-allele[as.character(loci_var),]
    freq2<-freq[as.character(loci_var),]
    allele_both<-rep(NA,length(loci_var))

    for(i in 1:length(loci_var)){
      var_index <- (mutation_collect[i,]>0)
      # identify the alleles at the variation locus
      allele_all <- allele2[i,]
      allele_var <- allele_all[var_index]
      allele_var <- allele_var[!is.na(allele_var)]

      # identify the corresponding freqencies of the alleles at the variation loci
      freq_all<-freq2[i,]
      freq_var<-freq_all[var_index]
      freq_var <- freq_var[!is.na(freq_var)]

      allele_all_var <- unlist(strsplit(allele_var,split="/"), use.names = F)
      freq_all_var <- as.numeric(unlist(strsplit(freq_var ,split="/"), use.names = F))

      # only includes alleles which have frequency within the [thre.lower, 1] interval
      allele_all_var <- allele_all_var[ freq_all_var >= thre.lower ]

      # combine with the reference allele
      if( sum(mutation_collect[i,] <= 1, na.rm = T) > 0 ){
        allele_all_var <- c( .mtRef[loci_var[i] ], allele_all_var)
      }
      allele_all_var <- allele_all_var[ !duplicated(allele_all_var) ]
      allele_all_var <- paste(allele_all_var,collapse = '/')
      allele_both[i] <- allele_all_var
    }

    head         <- c("mtID","ref_allele","allele_var","n_var","n_heter","n_homo","mut_allele",annot.select)

    #Open file for writing
    temp.file.name <- tempfile("var",fileext=c(".csv"))
    file.conn <- file( temp.file.name, "w")
    #file.conn <- file( paste0(path,study,"_annotation_variation.csv"), "w")
    write (head, file= file.conn ,sep=",",ncolumns=length(head), append=T)

    # this loop used to annote all the mutations at each heter loci with the scores chosen(annot.select)
    for (i in 1:length(loci_var)) {

      pos       <- loci_var[i]
      ref_allele <- .mtRef[pos]

      x    <- allele_both[i]
      x2   <- unlist(strsplit(x,split="/"), use.names = F)

      if(length(x2)>1){

        pos.not  <- which( !(x2 == ref_allele ))

        for (k in 1:length(pos.not)){
          point    <- x2[pos.not[k]]
          score = switch(point,
                         "A" = .AA[.AA$Pos==pos,annot.select],
                         "C" = .CC[.CC$Pos==pos,annot.select],
                         "G" = .GG[.GG$Pos==pos,annot.select],
                         "T" = .TT[.TT$Pos==pos,annot.select] )
        cat(paste(pos,ref_allele,x,n_mutation[i],heter_loci[loci==pos],homo_loci[loci==pos],point ,paste(score, collapse=","), sep=","),
              file=file.conn, fill=TRUE, append=TRUE)

        }
      } else if(length(x2)==1 & x2!=ref_allele){
        score = switch(x2,
                       "A" = .AA[.AA$Pos==pos,annot.select],
                       "C" = .CC[.CC$Pos==pos,annot.select],
                       "G" = .GG[.GG$Pos==pos,annot.select],
                       "T" = .TT[.TT$Pos==pos,annot.select] )
        cat(paste(pos,ref_allele,x,n_mutation[i],heter_loci[loci==pos],homo_loci[loci==pos],x2 ,paste(score, collapse=","), sep=","),
            file=file.conn, fill=TRUE, append=TRUE)
      }


    }

    close(file.conn)
    # Copy file to the permanent location:
    out.file.name <- paste0(path,study,"_annotation_variation.csv")
    file.copy(temp.file.name, out.file.name , overwrite = T)

  } else if(type=="heter"){

    allele2<-allele[as.character(loci_heter), ]
    freq2<-freq[as.character(loci_heter),]
    heter_collect<-mutation_collect[loci_var%in%loci_heter,]
    allele_heter<-rep(NA,length(loci_heter))

    for(i in 1:length(loci_heter)){
      heter_index<-(heter_collect[i,]==1)
      # identify the alleles at the variation locus
      allele_all<-allele2[i,]
      allele_var<-allele_all[heter_index]
      allele_var <- allele_var[!is.na(allele_var)]

      # identify the corresponding freqencies of the alleles at the variation loci
      freq_all<-freq2[i,]
      freq_var<-freq_all[heter_index]
      freq_var <- freq_var[!is.na(freq_var)]

      allele_all_var<-unlist(strsplit(allele_var,split="/"), use.names = F)
      freq_all_var<-as.numeric(unlist(strsplit(freq_var ,split="/"), use.names = F))
      # only includes alleles which have frequency within the [thre.lower, thre.upper] interval
      allele_all_var<-allele_all_var[freq_all_var>=thre.lower & freq_all_var<=thre.upper]
      # combine with the reference allele
      if(sum(heter_collect[i,]<=1, na.rm = T)>0){
        allele_all_var<-c(.mtRef[loci_heter[i]],allele_all_var)
      }
      allele_all_var<-allele_all_var[!duplicated(allele_all_var)]
      allele_all_var<-paste(allele_all_var,collapse = '/')
      allele_heter[i]<-allele_all_var
    }

    head         <- c("mtID","ref_allele","allele_heter","n_heter","mut_allele",annot.select)

    #Open file for writing
    temp.file.name <- tempfile("heter",fileext=c(".csv"))
    file.conn <- file( temp.file.name, "w")
    #file.conn <- file( paste0(path,study,"_annotation_heteroplasmy.csv"), "w")

    # output the annotation as a .csv file to the path users provided
    write (head, file= file.conn ,sep=",",ncolumns=length(head), append=T)

    # this loop used to annote all the heter mutations at each heter loci with the scores chosen(annot.select)
    for (i in 1:length(loci_heter)) {

      pos       <- loci_heter[i]
      ref_allele <- .mtRef[pos]

      x    <- allele_heter[i]
      x2   <- unlist(strsplit(x,split="/"), use.names = F)

      if(length(x2)>1){

        pos.not  <- which(!(x2==ref_allele))

        for (k in 1:length(pos.not)){
          point    <- x2[pos.not[k]]
          score = switch(point,
                         "A" = .AA[.AA$Pos==pos,annot.select],
                         "C" = .CC[.CC$Pos==pos,annot.select],
                         "G" = .GG[.GG$Pos==pos,annot.select],
                         "T" = .TT[.TT$Pos==pos,annot.select] )
          all   <- cbind(pos,ref_allele,x,heter_loci[loci==pos],point ,score)
          write.table(all, file= file.conn, sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
      }else if(length(x2)==1 & x2!=ref_allele){
        score = switch(x2,
                       "A" = .AA[.AA$Pos==pos,annot.select],
                       "C" = .CC[.CC$Pos==pos,annot.select],
                       "G" = .GG[.GG$Pos==pos,annot.select],
                       "T" = .TT[.TT$Pos==pos,annot.select] )
        all   <- cbind(pos,ref_allele,x,heter_loci[loci==pos],x2 ,score)
        write.table(all, file= file.conn, sep=",",row.names=F,col.names=F, quote=F, append=T)
      }
    }
    close(file.conn)
    # Copy file to the permanent location:
    out.file.name <- paste0(path,study,"_annotation_heteroplasmy.csv")
    file.copy(temp.file.name, out.file.name , overwrite = T)

  } else if(type=="homo"){

    allele2<-allele[as.character(loci_homo),]
    freq2<-freq[as.character(loci_homo),]
    homo_collect<-mutation_collect[loci_var%in%loci_homo,]
    allele_homo<-rep(NA,length(loci_homo))

    for(i in 1:length(loci_homo)){
      homo_index<-(homo_collect[i,]==2)
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
      if(sum(homo_collect[i,]<=1, na.rm = T)>0){
        allele_all_var<-c(.mtRef[loci_homo[i]],allele_all_var)
      }
      allele_all_var<-allele_all_var[!duplicated(allele_all_var)]
      allele_all_var<-paste(allele_all_var,collapse = '/')
      allele_homo[i]<-allele_all_var
    }

    head         <- c("mtID","ref_allele","allele_homo","n_homo","mut_allele",annot.select)

    #Open file for writing
    temp.file.name <- tempfile("homo",fileext=c(".csv"))
    file.conn <- file( temp.file.name, "w")
    #file.conn <- file( paste0(path,study,"_annotation_homoplasmy.csv"), "w")

    # output the annotation as a .csv file to the path users provided
    write (head, file= file.conn,sep=",",ncolumns=length(head), append=T)

    # this loop used to annote all the heter mutations at each heter loci with the scores chosen(annot.select)
    for (i in 1:length(loci_homo)) {

      pos       <- loci_homo[i]
      ref_allele <- .mtRef[pos]

      x    <- allele_homo[i]
      x2   <- unlist(strsplit(x,split="/"))

      if(length(x2)>1){

        pos.not  <- which(!(x2==ref_allele))
        for (k in 1:length(pos.not)){
          point    <- x2[pos.not[k]]
          score = switch(point,
                         "A" = .AA[.AA$Pos==pos,annot.select],
                         "C" = .CC[.CC$Pos==pos,annot.select],
                         "G" = .GG[.GG$Pos==pos,annot.select],
                         "T" = .TT[.TT$Pos==pos,annot.select] )
          all   <- cbind(pos,ref_allele,x,homo_loci[loci==pos],point ,score)
          write.table(all, file= file.conn,sep=",",row.names=F,col.names=F, quote=F, append=T)
        }
      } else if(length(x2)==1 & x2!=ref_allele){
        score = switch(x2,
                       "A" = .AA[.AA$Pos==pos,annot.select],
                       "C" = .CC[.CC$Pos==pos,annot.select],
                       "G" = .GG[.GG$Pos==pos,annot.select],
                       "T" = .TT[.TT$Pos==pos,annot.select] )
        all   <- cbind(pos,ref_allele,x,homo_loci[loci==pos],x2 ,score)
        write.table(all, file= file.conn,sep=",",row.names=F,col.names=F, quote=F, append=T)
      }
    }
    close(file.conn)
    # Copy file to the permanent location:
    out.file.name <- paste0(path,study,"_annotation_homoplasmy.csv")
    file.copy(temp.file.name, out.file.name , overwrite = T)
  }
  return(mt_summary_obj)
}










