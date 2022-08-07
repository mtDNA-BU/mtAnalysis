#' mtSummary function
#'
#' @description Identification of mtDNA variations, output summary statistics,
#' annotation of heteroplasmic and/or homoplasmic variations.
#'
#' @param aaf a numeric matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' It contains subject ID as the column names,
#' and the AAFs of all 16569 mtDNA loci for each subject.
#' It is generated from mtAAF function.
#' @param allele a character matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains the alleles of each subject at each locus.
#' The matrix must contain subject ID as the column names.
#' "/" is used to delimited different allele calls in a locus.
#' @param freq a character matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains the allele fractions of the corresponding allele matrix.
#' The matrix must contain subject ID as the column names.
#' "/" is used to delimited the allele fractions.
#' @param coverage a numeric matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains the reads coverage of the 16569 mtDNA loci for each
#' subject. The matrix must contain the subject ID as the column names.
#' @param coverage.qc a number(default is 250) of threshold for the coverage.
#' If the coverage<coverage.qc, the allele call at that locus of the subject
#' will not be used.
#' @param thre.lower a number(default is 0.03) of lower bound of the threshold
#' defining heteroplasmic
#' and homoplasmic variations
#' @param thre.upper a number(default is 0.97) of upper bound of the threshold
#' defining heteroplasmic
#' and homoplasmic variations
#' @param loci one of:
#' 1. a vector(default is c(1:16569)) of mitochondrial DNA  loci to specify
#' which loci should be used to identify the variations and annotate,
#' 2. a character string for the regions
#' (e.g. "coding" , "tRNA", "RNR1" , "RNR2",...)
#' @param type a character of indicator choosing to output annotation
#' to all variations,
#' heteroplasmic variations, or homoplasmic variations. “both” returns
#' annotation to all
#' variations (default), "heter" returns annotation to heteroplasmic
#' variations and "homo"
#' returns annotation to homoplasmic variations.
#' @param coverSummary logical(default is True). A user can specify to output
#' summary of mean coverage at each mtDNA loci (across all participants) and
#' for each individual (across all mtDNA loci).
#' @param varHist logical(default is True). A user can specify to output
#' histograms to visualize the
#' heteroplasmic and homoplasmic burden across participants and mtDNA loci.
#' @param annot.select A character vector of variation position, alternative
#' allele, corresponding gene and types of annotation scores to output based on
#' user's choice. The available choices are "Pos", "ref", "Gene", "TypeMutation",
#' "MissensMutation", "CodonPosition", "ProteinDomain", "mFOLD_dG",
#' "mFOLD_Initial", "mFOLD_rCRS.DG", "mFOLD_rCRS.Initial",
#' "mFOLD_AnticodonAminoAcidChange", "mFOLD_Location", "PolyPhen2",
#' "PolyPhen2_score", "SIFT", "SIFT_score", "PROVEAN", "PROVEAN_score",
#' "MutationAssessor", "MutationAssessor_score", "CADD", "CADD_score",
#' "CADD_phred_score", "PANTHER", "PANTHER_score", "PhD_SNP", "PhD_SNP_score",
#' "SNAP", "SNAP_score", "MutationTaster", "MutationTaster_score", "dbSNP_150_id",
#' "APOGEE_boost_consensus", "APOGEE_boost_mean_prob", "APOGEE_boost_mean"
#' @param path the path of the output annotation file. If not provided,
#' the annotation file will
#' output to the current working directory
#' @param study A string of study names. Default is Study.
#' @param anno A logical value (default is False) indicating whether output the
#' annotation results in a .csv file
#'
#' @return 1. Summary frequency: descriptive statistics of sequencing coverage
#' across individual
#' and across mtDNA loci; the total number of mtDNA loci with variations and
#' the number of
#' heteroplasmic/homoplasmic mtDNA loci in the study sample;
#' the min/Q1/mean/median/Q3/largest number
#' of homoplasmy/heteroplasmy carried by an individual.
#' 2. Plots: a scatter plot of the median coverage
#' across mtDNA loci; histograms to visualize the heteroplasmic and
#' homoplasmic burden across
#' participants and mtDNA loci based on user’s choice.
#' 3.	Summary annotation for all of the
#' heteroplasmic and/or homoplasmic variations observed in
#' the study data based on user's choice.
#' @import graphics grDevices utils
#' @export
#' @examples
#'
#'\dontrun{
#' ## Read input data
#' allele_file <- "allele.csv"
#' freq_file   <- "freq.csv"
#'
#' allele <- as.matrix( read.csv(file = allele_file, sep = ",") )
#' freq   <- as.matrix( read.csv(file = freq_file, sep = ",") )
#'
#' aaf =  mtAAF ( allele, freq)
#'
#' ## Summary for the aaf object
#' mtSummary(aaf, allele, freq, coverage)
#'}

mtSummary<-function(aaf, allele, freq, coverage,
                    coverage.qc=250,
                    thre.lower=0.03, thre.upper=0.97,
                    loci=c(1 : .mtLength),
                    type="both",
                    coverSummary=T, varHist=T,
                    annot.select=c("Pos", "ref", "Gene", "TypeMutation",
                                   "MissensMutation", "CodonPosition",
                                   "ProteinDomain", "dbSNP_150_id", "PolyPhen2",
                                   "PolyPhen2_score",
                                   "SIFT", "SIFT_score", "CADD", "CADD_score",
                                   "CADD_phred_score"),
                    path="./",
                    study="Study",
                    anno=T) {

    if(!all(is.character(allele), is.character(freq) ,
            is.numeric(aaf), is.numeric(coverage)))
        stop("the allele and freq shoud be character,
             aaf and coverage should be numeric")


    if(!all(sapply(list(dim(aaf), dim(allele), dim(freq)),
                   function(x) x == dim(coverage))))
        stop("the coverage, allele, frequency and
             AAF should have same dimension")


    if(dim(aaf)[1] != .mtLength)
        stop("the coverage, allele, frequency and aaf should have
             16569 loci (rows)")


    if(!all(setequal(colnames(allele), colnames(freq)),
            setequal(colnames(allele), colnames(aaf)),
            setequal(colnames(allele), colnames(coverage))))
        stop("the coverage, allele, frequency and aaf should have same
             subject IDs (column names)")

    if(!(is.logical(coverSummary) & length(coverSummary)==1)){
        stop("coverSummary must be logical value")
    }

    if(!(is.logical(varHist) & length(varHist)==1)){
        stop("varHist must be logical value")
    }

    if(!(is.logical(anno) & length(anno)==1)){
        stop("anno must be logical value")
    }


    if(!all(mode(loci) %in% c("numeric","character"), is.vector(loci)))
        stop("loci must be numeric vector or character")


    if(is.character(loci) & !all(loci %in% c("coding", "tRNA", "RNR1", "RNR2") &
         length(loci) == 1))
        stop("loci must be one of coding, tRNA, RNR1, RNR2 if it is character")


    if(is.numeric(loci) & !all(loci %in% (1:.mtLength)))
        stop("loci should be a subset of 1:16569")


    if( ! dir.exists(path) )
        stop( paste(" Output path", path, "does not exist.") )

    if(!all(is.numeric(coverage.qc), length(coverage.qc)==1))
        stop("coverage.qc must be numeric of length 1")


    if(!all(is.numeric(c(thre.lower, thre.upper)),
            length(thre.lower) == 1,
            length(thre.upper) == 1))
        stop("thre.lower and thre.upper must be numeric of length 1")


    if(thre.lower >= thre.upper)
        stop("thre.lower must be smaller than thre.upper")


    if(!all(length(type) == 1 , type %in% c("both","heter","homo")))
        stop("type must be one of both, heter and homo")


    ## Check if annot.select belongs to the colnames of annotation files
    if(!all(annot.select %in% names(.AA), is.vector(annot.select)))
        stop("annot.select must be a subset of variables of annotation files")


    subjectID <- colnames(allele)
    subjectID <- sort(subjectID)
    allele <- allele[ , subjectID]
    freq <- freq[ , subjectID]
    coverage <- coverage[ , subjectID]
    aaf <- aaf[ , subjectID]

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

    ## categorize by AAF: AAF < thre.lower means no variation(coded as 0)
    ## thre.lower <= AAF <= thre.upper means heteroplasmic variation(coded as 1)
    ## thre.upper<AAF means homoplastic variation(coded as 2)
    aaf_cat <- matrix(0, nrow=nrow(aaf), ncol=ncol(aaf))
    aaf_cat[ aaf >= thre.lower & aaf <= thre.upper ] <- 1
    aaf_cat[ aaf > thre.upper ] <- 2
    rownames(aaf_cat) <- rownames(aaf)
    colnames(aaf_cat) <- colnames(aaf)


    aaf_cat[ coverage < coverage.qc ] <- NA

    ## The reads of these loci are not reliable, so these loci should be removed
    loci_removed <- .mtLociNUMT   # c(301,302,310,316,3107,16182)
    loci_removed <- as.character(loci_removed)
    loci_left <- setdiff(rownames(aaf_cat), loci_removed)
    aaf_cat <- aaf_cat[loci_left,]

    ## only include mutations loci: at least one subject has mutation at that locus
    mutation_collect <- aaf_cat

    n_mutation <- rowSums(aaf_cat > 0, na.rm=T)

    ## only include the mutation(heter/homo) loci
    mutation_collect <- mutation_collect[n_mutation > 0,]

    ## loci of variations(heter/homo)
    loci_var <- as.numeric(rownames(mutation_collect))
    mt_summary_obj$loci_var <- loci_var
    n_mutation <- n_mutation[n_mutation > 0]

    ## output IDs of individuals
    mt_summary_obj$IDs <- subjectID

    ## output summary statistics of hereoplasmic variations
    ## calculate the number of heteroplasmic mutations for each subject
    heter_burden  <- colSums(mutation_collect == 1, na.rm=T)
    mt_summary_obj$heter_burden <- heter_burden
    mt_summary_obj$heter_burden_sum <- summary(heter_burden)

    ## calculate the number of heteroplasmic mutations for each locus
    heter_loci  <- rowSums(aaf_cat == 1, na.rm=T)
    mt_summary_obj$heter_loci <- heter_loci
    mt_summary_obj$heter_loci_sum <- summary(heter_loci)

    ## loci with heteroplasmic variations
    loci_heter <- as.numeric(names(heter_loci))[heter_loci > 0]
    mt_summary_obj$loci_heter <- loci_heter

    ## calculate the total number of heteroplasmic variations
    mt_summary_obj$heter_total <- sum(heter_burden)

    ## output summary statistics of homoplasmic variations
    ## calculate the number of homoplasmic variations for each subject
    homo_burden  <- colSums(mutation_collect == 2, na.rm=T)
    mt_summary_obj$homo_burden <- homo_burden
    mt_summary_obj$homo_burden_sum <- summary(homo_burden)

    ## calculate the number of homoplasmic mutations for each locus
    homo_loci  <- rowSums(aaf_cat == 2, na.rm=T)
    mt_summary_obj$homo_loci <- homo_loci
    mt_summary_obj$homo_loci_sum <- summary(homo_loci)

    ## loci with homoplasmic variations
    loci_homo <- as.numeric(names(homo_loci))[homo_loci > 0]
    mt_summary_obj$loci_homo <- loci_homo

    ## calculate the total number of homoplasmic variations
    mt_summary_obj$homo_total <- sum(homo_burden)

    ## generate histograms of heter and homo burden of across subjects and across mtDNA loci
    if(varHist){

        ## Open an output pdf file
        pdf(file = paste0(path, "/mtHistograms.pdf"))

        ## histogram of heteroplasmic burden across subjects
        h <- hist(heter_burden,
                  breaks = 100,
                  xlab ="Heteroplasmic burden score" ,
                  main ="Histogram of heteroplasmic burden score" )

        h <- hist(heter_loci,
                  breaks = 100,
                  xlab ="Heteroplasmic variations of mtDNA loci" ,
                  main ="Histogram of heteroplasmic variations of mtDNA loci" )

        h <- hist(homo_burden,
                  breaks = 100,
                  xlab ="Homoplasmic burden score" ,
                  main ="Histogram of homoplasmic burden score" )

        h <- hist(homo_loci,
                  breaks = 100,
                  xlab ="Homoplasmic variations of mtDNA loci" ,
                  main ="Histogram of homoplasmic variations of mtDNA loci" )

        dev.off()
    }

    if(anno){
    ## annotation for all heter and/or homo variations at each
    ## mutation loci based on user's choice
    ## if type=="both", annotate both heter/homo variations
    if(type == "both"){

        allele2 <- allele[as.character(loci_var),]
        freq2 <- freq[as.character(loci_var),]
        allele_both <- rep(NA,length(loci_var))

        for(i in seq_along(loci_var)) {
            var_index <- (mutation_collect[i,] > 0)

            ## identify the alleles at the variation locus
            allele_all <- allele2[i,]
            allele_var <- allele_all[var_index]
            allele_var <- allele_var[ !is.na(allele_var)]

            ## identify the corresponding frequencies of
            ## the alleles at the variation loci
            freq_all<-freq2[i,]
            freq_var<-freq_all[var_index]
            freq_var <- freq_var[ !is.na(freq_var)]

            allele_all_var <- unlist(strsplit(allele_var,split="/"),
                                     use.names = F)
            freq_all_var <- as.numeric(unlist(strsplit(freq_var ,
                                                       split="/"),
                                              use.names = F))

            ## only includes alleles which have frequency
            ## within the [thre.lower, 1] interval
            allele_all_var <- allele_all_var[ freq_all_var >= thre.lower ]

            ## combine with the reference allele
            if( sum(mutation_collect[i,] == 0, na.rm = T) > 0 ){
                allele_all_var <- c( .mtRef[loci_var[i] ], allele_all_var)
            }
            allele_all_var <- allele_all_var[ !duplicated(allele_all_var) ]
            allele_all_var <- paste(allele_all_var,collapse = '/')
            allele_both[i] <- allele_all_var
        }

        head <- c("mtID", "ref_allele", "allele_all", "n_var",
                  "n_heter", "n_homo", "mut_allele", annot.select)

        ## Open file for writing
        temp.file.name <- tempfile("var",fileext=c(".csv"))
        file.conn <- file( temp.file.name, "w")
        write (head, file= file.conn ,sep=",",ncolumns=length(head), append=T)

        ## this loop used to annotate all the mutations at each
        ## heter loci with the scores chosen(annot.select)
        for (i in seq_along(loci_var)) {

            pos <- loci_var[i]
            ref_allele <- .mtRef[pos]

            x    <- allele_both[i]
            x2   <- unlist(strsplit(x,split="/"), use.names = F)

            if(length(x2) > 1){

                pos.not  <- which( !(x2 == ref_allele ))
                if(length(pos.not)==1){
                    point    <- x2[pos.not]
                    score = switch(point,
                                   "A" = .AA[.AA$Pos == pos,annot.select],
                                   "C" = .CC[.CC$Pos == pos,annot.select],
                                   "G" = .GG[.GG$Pos == pos,annot.select],
                                   "T" = .TT[.TT$Pos == pos,annot.select] )
                    if(point=="D" | point=="I"){
                        score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                        NA_len <- length(annot.select)-length(score)
                        score <- c(score, rep(NA, NA_len))
                    }
                    cat(paste(pos,
                              ref_allele,
                              x,
                              n_mutation[i],
                              heter_loci[as.character(pos)],
                              homo_loci[as.character(pos)],
                              point,
                              paste(score, collapse=","), sep=","),
                        file=file.conn, fill=TRUE, append=TRUE)


                }else if(length(pos.not)>1){
                    var_index <- (mutation_collect[i,] > 0)

                    ## identify the alleles at the variation locus
                    allele_all <- allele2[i,]
                    allele_var <- allele_all[var_index]
                    allele_var <- allele_var[ !is.na(allele_var)]

                    ## identify the corresponding frequencies of
                    ## the alleles at the variation loci
                    freq_all<-freq2[i,]
                    freq_var<-freq_all[var_index]
                    freq_var <- freq_var[ !is.na(freq_var)]

                    allele_all_var <- unlist(strsplit(allele_var,split="/"),
                                             use.names = F)
                    freq_all_var <- as.numeric(unlist(strsplit(freq_var ,
                                                               split="/"),
                                                      use.names = F))

                    ## only includes alleles which have frequency
                    ## within the [thre.lower, 1] interval
                    allele_all_var_heter <- allele_all_var[ freq_all_var >= thre.lower & freq_all_var<=thre.upper]
                    allele_all_var_homo <- allele_all_var[ freq_all_var > thre.upper ]

                for (k in seq_along(pos.not)){
                    point    <- x2[pos.not[k]]
                    score = switch(point,
                                   "A" = .AA[.AA$Pos == pos,annot.select],
                                   "C" = .CC[.CC$Pos == pos,annot.select],
                                   "G" = .GG[.GG$Pos == pos,annot.select],
                                   "T" = .TT[.TT$Pos == pos,annot.select])
                    if(point=="D" | point=="I"){
                        score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                        NA_len <- length(annot.select)-length(score)
                        score <- c(score, rep(NA, NA_len))
                    }
                    n_heter1 <- sum(point == allele_all_var_heter)
                    n_homo1 <- sum(point == allele_all_var_homo)
                    n_var1 <- n_heter1 + n_homo1
                    cat(paste(pos,
                              ref_allele,
                              x,
                              n_var1,
                              n_heter1,
                              n_homo1,
                              point,
                              paste(score, collapse=","), sep=","),
                        file=file.conn, fill=TRUE, append=TRUE)

                }
                }
            } else if(length(x2) == 1 & x2 != ref_allele){
                score = switch(x2,
                               "A" = .AA[.AA$Pos == pos,annot.select],
                               "C" = .CC[.CC$Pos == pos,annot.select],
                               "G" = .GG[.GG$Pos == pos,annot.select],
                               "T" = .TT[.TT$Pos == pos,annot.select])
                if(point=="D" | point=="I"){
                    score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                    NA_len <- length(annot.select)-length(score)
                    score <- c(score, rep(NA, NA_len))
                }
                cat(paste(pos,
                          ref_allele,
                          x,
                          n_mutation[i],
                          heter_loci[as.character(pos)],
                          homo_loci[as.character(pos)],
                          x2,
                          paste(score, collapse=","), sep=","),
                    file=file.conn, fill=TRUE, append=TRUE)
            }

        }

        close(file.conn)


        ## Copy file to the permanent location:
        out.file.name <- paste0(path, study, "_annotation_variation.csv")
        file.copy(temp.file.name, out.file.name , overwrite = T)

    } else if(type == "heter"){

        allele2 <- allele[as.character(loci_heter), ]
        freq2 <- freq[as.character(loci_heter),]
        heter_collect <- mutation_collect[loci_var %in% loci_heter,]
        allele_heter <- rep(NA, length(loci_heter))

        for(i in seq_along(loci_heter)){
            heter_index <- (heter_collect[i,] == 1)

            ## identify the alleles at the variation locus
            allele_all <- allele2[i,]
            allele_var <- allele_all[heter_index]
            allele_var <- allele_var[ !is.na(allele_var)]

            ## identify the corresponding frequencies of the alleles at the variation loci
            freq_all <- freq2[i,]
            freq_var <- freq_all[heter_index]
            freq_var <- freq_var[ !is.na(freq_var)]

            allele_all_var <- unlist(strsplit(allele_var,split="/"),
                                     use.names = F)

            freq_all_var <- as.numeric(unlist(strsplit(freq_var ,split="/"),
                                            use.names = F))

            ## only includes alleles which have frequency within the [thre.lower, thre.upper] interval
            allele_all_var <- allele_all_var[freq_all_var >= thre.lower &
                                             freq_all_var <= thre.upper]

            ## combine with the reference allele
            if(sum(heter_collect[i,] == 0, na.rm = T) > 0){
                allele_all_var <- c(.mtRef[loci_heter[i]], allele_all_var)
            }

            allele_all_var <- allele_all_var[!duplicated(allele_all_var)]
            allele_all_var <- paste(allele_all_var,collapse = '/')
            allele_heter[i] <- allele_all_var
        }

        head <- c("mtID", "ref_allele", "allele_all", "n_heter", "mut_allele",
                  annot.select)

        ## Open file for writing
        temp.file.name <- tempfile("heter", fileext=c(".csv"))
        file.conn <- file( temp.file.name, "w")
        write (head, file=file.conn, sep=",", ncolumns=length(head), append=T)

        ## this loop used to annote all the heter mutations at each heter loci
        ## with the scores chosen(annot.select)
        for (i in seq_along(loci_heter)) {

            pos <- loci_heter[i]
            ref_allele <- .mtRef[pos]

            x    <- allele_heter[i]
            x2   <- unlist(strsplit(x,split="/"), use.names = F)

            if(length(x2) > 1){

                pos.not <- which(!(x2 == ref_allele))
                if(length(pos.not)==1){
                    point    <- x2[pos.not]
                    score = switch(point,
                                   "A" = .AA[.AA$Pos == pos,annot.select],
                                   "C" = .CC[.CC$Pos == pos,annot.select],
                                   "G" = .GG[.GG$Pos == pos,annot.select],
                                   "T" = .TT[.TT$Pos == pos,annot.select])
                    if(point=="D" | point=="I"){
                        score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                        NA_len <- length(annot.select)-length(score)
                        score <- c(score, rep(NA, NA_len))
                    }
                    cat(paste(pos,
                              ref_allele,
                              x,
                              heter_loci[as.character(pos)],
                              point,
                              paste(score, collapse=","), sep=","),
                        file=file.conn, fill=TRUE, append=TRUE)

                }else if(length(pos.not)>1){
                    heter_index <- (heter_collect[i,] == 1)

                    ## identify the alleles at the variation locus
                    allele_all <- allele2[i,]
                    allele_var <- allele_all[heter_index]
                    allele_var <- allele_var[ !is.na(allele_var)]

                    ## identify the corresponding frequencies of the alleles at the variation loci
                    freq_all <- freq2[i,]
                    freq_var <- freq_all[heter_index]
                    freq_var <- freq_var[ !is.na(freq_var)]

                    allele_all_var <- unlist(strsplit(allele_var,split="/"),
                                             use.names = F)

                    freq_all_var <- as.numeric(unlist(strsplit(freq_var ,split="/"),
                                                      use.names = F))

                    ## only includes alleles which have frequency
                    ## within the [thre.lower, thre.upper] interval
                    allele_all_var <- allele_all_var[freq_all_var >= thre.lower &
                                                         freq_all_var <= thre.upper]

                for (k in seq_along(pos.not)){
                    point    <- x2[pos.not[k]]
                    score = switch(point,
                                   "A" = .AA[.AA$Pos == pos,annot.select],
                                   "C" = .CC[.CC$Pos == pos,annot.select],
                                   "G" = .GG[.GG$Pos == pos,annot.select],
                                   "T" = .TT[.TT$Pos == pos,annot.select])
                    if(point=="D" | point=="I"){
                        score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                        NA_len <- length(annot.select)-length(score)
                        score <- c(score, rep(NA, NA_len))
                    }
                    cat(paste(pos,
                              ref_allele,
                              x,
                              sum(point == allele_all_var),
                              point,
                              paste(score, collapse=","), sep=","),
                        file=file.conn, fill=TRUE, append=TRUE)

                }
                }
            }else if(length(x2) == 1 & x2 != ref_allele){
                score = switch(x2,
                               "A" = .AA[.AA$Pos == pos,annot.select],
                               "C" = .CC[.CC$Pos == pos,annot.select],
                               "G" = .GG[.GG$Pos == pos,annot.select],
                               "T" = .TT[.TT$Pos == pos,annot.select])
                if(point=="D" | point=="I"){
                    score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                    NA_len <- length(annot.select)-length(score)
                    score <- c(score, rep(NA, NA_len))
                }
                cat(paste(pos,
                          ref_allele,
                          x,
                          heter_loci[as.character(pos)],
                          x2,
                          paste(score, collapse=","), sep=","),
                    file=file.conn, fill=TRUE, append=TRUE)

            }
        }

        close(file.conn)

        ## Copy file to the permanent location:
        out.file.name <- paste0(path, study, "_annotation_heteroplasmy.csv")
        file.copy(temp.file.name, out.file.name , overwrite = T)

    } else if(type=="homo"){

        allele2 <- allele[as.character(loci_homo),]
        freq2 <- freq[as.character(loci_homo),]
        homo_collect <- mutation_collect[loci_var %in% loci_homo,]
        allele_homo <- rep(NA,length(loci_homo))

        for(i in seq_along(loci_homo)){
            homo_index <- (homo_collect[i,]==2)
            ## identify the alleles at the variation locus
            allele_all <- allele2[i,]
            allele_var <- allele_all[homo_index]
            allele_var <- allele_var[!is.na(allele_var)]

            ## identify the corresponding freqencies of the alleles
            ## at the variation loci
            freq_all <- freq2[i,]
            freq_var <- freq_all[homo_index]
            freq_var <- freq_var[!is.na(freq_var)]

            allele_all_var <- unlist(strsplit(allele_var, split="/"))
            freq_all_var <- as.numeric(unlist(strsplit(freq_var, split="/")))

            ## only includes alleles which have frequency within
            ##the [thre.lower, thre.upper] interval
            allele_all_var<-allele_all_var[freq_all_var > thre.upper]

            ## combine with the reference allele
            if(sum(homo_collect[i,] == 0, na.rm = T) > 0){
                allele_all_var <- c(.mtRef[loci_homo[i]], allele_all_var)
            }
            allele_all_var <- allele_all_var[!duplicated(allele_all_var)]
            allele_all_var <- paste(allele_all_var,collapse = '/')
            allele_homo[i] <- allele_all_var
        }

        head <- c("mtID", "ref_allele", "allele_all", "n_homo", "mut_allele",
                  annot.select)

        #Open file for writing
        temp.file.name <- tempfile("homo",fileext=c(".csv"))
        file.conn <- file( temp.file.name, "w")
        write (head, file= file.conn, sep=",", ncolumns=length(head), append=T)

        ## this loop used to annote all the homo mutations at each
        ## homo loci with the scores chosen(annot.select)
        for (i in seq_along(loci_homo)) {

            pos <- loci_homo[i]
            ref_allele <- .mtRef[pos]

            x    <- allele_homo[i]
            x2   <- unlist(strsplit(x,split="/"))

            if(length(x2)>1){

                pos.not  <- which(!(x2==ref_allele))
                if(length(pos.not)==1){
                    point    <- x2[pos.not]
                    score = switch(point,
                                   "A" = .AA[.AA$Pos == pos,annot.select],
                                   "C" = .CC[.CC$Pos == pos,annot.select],
                                   "G" = .GG[.GG$Pos == pos,annot.select],
                                   "T" = .TT[.TT$Pos == pos,annot.select])
                    if(point=="D" | point=="I"){
                        score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                        NA_len <- length(annot.select)-length(score)
                        score <- c(score, rep(NA, NA_len))
                    }
                    cat(paste(pos,
                              ref_allele,
                              x,
                              homo_loci[as.character(pos)],
                              point,
                              paste(score, collapse=","), sep=","),
                        file=file.conn, fill=TRUE, append=TRUE)

                }else if(length(pos.not)>1){
                    homo_index <- (homo_collect[i,]==2)
                    ## identify the alleles at the variation locus
                    allele_all <- allele2[i,]
                    allele_var <- allele_all[homo_index]
                    allele_var <- allele_var[!is.na(allele_var)]

                    ## identify the corresponding freqencies of the alleles
                    ## at the variation loci
                    freq_all <- freq2[i,]
                    freq_var <- freq_all[homo_index]
                    freq_var <- freq_var[!is.na(freq_var)]

                    allele_all_var <- unlist(strsplit(allele_var, split="/"))
                    freq_all_var <- as.numeric(unlist(strsplit(freq_var, split="/")))

                    ## only includes alleles which have frequency within
                    ## the (thre.upper, 1] interval
                    allele_all_var<-allele_all_var[freq_all_var > thre.upper]

                    for (k in seq_along(pos.not)){
                        point    <- x2[pos.not[k]]
                        score = switch(point,
                                       "A" = .AA[.AA$Pos == pos,annot.select],
                                       "C" = .CC[.CC$Pos == pos,annot.select],
                                       "G" = .GG[.GG$Pos == pos,annot.select],
                                       "T" = .TT[.TT$Pos == pos,annot.select])
                        if(point=="D" | point=="I"){
                            score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                            NA_len <- length(annot.select)-length(score)
                            score <- c(score, rep(NA, NA_len))
                        }
                        cat(paste(pos,
                                  ref_allele,
                                  x,
                                  sum(point==allele_all_var),
                                  point,
                                  paste(score, collapse=","), sep=","),
                            file=file.conn, fill=TRUE, append=TRUE)
                    }
                }
            } else if(length(x2)==1 & x2 != ref_allele){
                score = switch(x2,
                               "A" = .AA[.AA$Pos == pos,annot.select],
                               "C" = .CC[.CC$Pos == pos,annot.select],
                               "G" = .GG[.GG$Pos == pos,annot.select],
                               "T" = .TT[.TT$Pos == pos,annot.select])
                if(point=="D" | point=="I"){
                    score <- .AA[.AA$Pos == pos,c("Pos", "ref", "Gene")]
                    NA_len <- length(annot.select)-length(score)
                    score <- c(score, rep(NA, NA_len))
                }
                cat(paste(pos,
                          ref_allele,
                          x,
                          homo_loci[as.character(pos)],
                          x2,
                          paste(score, collapse=","), sep=","),
                    file=file.conn, fill=TRUE, append=TRUE)
            }
        }
        close(file.conn)

        ## Copy file to the permanent location:
        out.file.name <- paste0(path,study, "_annotation_homoplasmy.csv")
        file.copy(temp.file.name, out.file.name , overwrite = T)
    }
    }

    return(mt_summary_obj)
}










