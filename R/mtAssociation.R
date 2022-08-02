# function to run association analysis of heteroplasmies





#' mtAssociation function
#'
#' @description Association analysis of heteroplasmies with disease traits by
#' four common methods of burden test, SKAT, SKAT-O and ACAT-O. ACAT-O is used
#' to combine results of burden test and SKAT.
#'
#' @param aaf a numeric matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' It contains subject ID as the column names,
#' and the AAFs of all 16569 mtDNA loci for each subject.
#' It is generated from mtAAF function.
#' @param coverage a numeric matrix (16569 x N) provided by the user.
#' Rows correspond to loci and columns correspond to subjects.
#' This matrix contains the reads coverage of the 16569 mtDNA loci for each
#' subject. The matrix must contain the subject ID as the column names.
#' @param coverage.qc a number(default is 250) of threshold for the coverage.
#' If the coverage<coverage.qc, the allele call at that locus of the subject
#' will not be used.
#' @param family a character string to specify data type of outcome: "gaussian"
#' for continuous data and "binomial" for binary data.
#' @param pheno a data frame containing data of covariates and outcome. Rows are
#' subjects and columns are variables. It must have the same sample size as the
#' aaf matrix. It must contain the subject ID as the row names.
#' @param trait a character string of outcome name to be analyzed.
#' @param covars a character vector indicating the names of covariates.
#' @param G_coding Coding definition of heteroplasmy. "heter" for coding
#' definition 1, and "heter2" for coding definition 2. See details.
#' @param rho_skatO A sequence of values that specify combinations of SKAT and a
#' burden test to be considered in SKAT-O (default is c(0,1)).
#' @param heter_scale logical(default is False). A user can specify to
#' standardize the genetic dosage of heteroplasmies of each variant.
#' @param region a numeric vector to specify the gene region to be analyzed
#' (default is c(1:16569)).
#' @param thre_lower a number (default is 0.03) of lower bound of the threshold
#' defining heteroplasmic mutations
#' @param thre_upper a number (default is 0.97) of upper bound of the threshold
#' defining heteroplasmic mutations
#' @param maf_max upper bound of frequency of heteroplasmies at population level
#' to be included
#' @param kins Kinship matrix for related subjects. Only support for continuous
#' outcomes currently.
#'
#' @details By coding definition 1 (G_coding="heter"), the genetic dosage of
#' heteroplasmy is 1 if the corresponding alternative allele fraction (aaf)
#' satisfies thre_lower≤aaf≤thre_upper, and 0 otherwise. By coding definition 2
#' (G_coding="heter2"), the genetic dosage of heteroplasmy is the aaf itself if
#' the corresponding alternative allele fraction (aaf) satisfies
#' thre_lower≤aaf≤thre_upper, and 0 otherwise.
#'
#' @return A data frame containing the results of association analysis of
#' heteroplasmies by the four methods.
#' @export
#'
#' @examples

mtAssociation<-function(aaf, coverage, coverage.qc=250, family, pheno,
                     trait, covars, G_coding="heter", rho_skatO=c(0,1), heter_scale=F,
                     region=c(1:16569), thre_lower=0.03, thre_upper=0.97, maf_max=0.01,
                     kins=NULL){

    ###input type check ###
    if (length(coverage.qc)!=1) stop("Argument coverage.qc must be of length one")
    if (length(family)!=1) stop("Argument family must be of length one")
    if (length(G_coding)!=1) stop("Argument G_coding must be of length one")
    if (length(heter_scale)!=1) stop("Argument heter_scale must be of length one")
    if (length(thre_lower)!=1) stop("Argument thre_lower must be of length one")
    if (length(thre_upper)!=1) stop("Argument thre_upper must be of length one")
    if (length(maf_max)!=1) stop("Argument maf_max must be of length one")

    if(!is.numeric(aaf))
        stop("Argument aaf is not a numeric matrix")
    if(!is.numeric(coverage))
        stop("Argument coverage is not a numeric matrix")
    if(is.null(colnames(aaf)))
        stop("Argument aaf must contain the subject ID as the column names")
    if(is.null(colnames(coverage)))
        stop("Argument coverage must contain the subject ID as the column names")

    if(!is.numeric(coverage.qc))
        stop("Argument coverage.qc is not numeric")
    if(coverage.qc<0)
        stop("Argument coverage.qc must be positive")

    if(!family %in% c("gaussian","binomial"))
        stop('Argument family must be a string input of "gaussian" or "binomial"')

    if(!G_coding %in% c("heter","heter2"))
        stop('Argument family must be a string input of "heter" or "heter2"')

    if(!is.numeric(region))
        stop("Argument region is not a numeric vector")
    if(sum(region!=round(region))>0)
        stop("Argument region must have integer entries")
    if(!is.numeric(thre_lower))
        stop("Argument thre_lower is not numeric")
    if(!is.numeric(thre_upper))
        stop("Argument thre_upper is not numeric")
    if(!is.numeric(maf_max))
        stop("Argument maf_max is not numeric")

    if(!is.numeric(rho_skatO))
        stop("Argument rho_skatO is not numeric")
    if(sum(duplicated(rho_skatO))!=0)
        stop("Argument rho_skatO contains duplicates")

    if(!is.logical(heter_scale))
        stop("Argument heter_scale is not logical")

    if (dim(pheno)[1] != dim(aaf)[2])
        stop("dataframe pheno must have the same sample size as the aaf matrix")
    if (sum(!colnames(aaf) %in% row.names(pheno)) >0)
        stop("pheno row names must have all subject IDs present in aaf column names")
    if(!is.character(colnames(pheno)))
        stop("Argument pheno must contain variabe names as column names")

    if(!is.character(covars))
        stop("Argument covars is not a character vector")
    if(!is.character(trait) || length(trait)!=1)
        stop("Argument trait is not a single character string")
    if(!covars %in% names(pheno))
        stop("Argument covars must match the covariate names in dataframe pheno")
    if(!trait %in% names(pheno))
        stop("Argument trait must match the outcome variable name in dataframe pheno")

    if(!is.null(kins)){
        if (sum(!colnames(aaf) %in% colnames(kins)) >0)
            stop("matrix kins must have all subject IDs present in aaf column names")
    }

    ### input range check ###

    if (min(region)<1 || max(region)>16569)
        stop("Argument region out of range 1:16569")

    if (thre_lower < 0 || thre_lower > 1)
        stop("Argument thre_lower out of range")
    if (thre_upper < 0 || thre_upper > 1)
        stop("Argument thre_upper out of range")
    if (thre_lower >= thre_upper)
        stop("Argument thre_lower must be smaller than thre_upper")

    if (maf_max < 0 || maf_max > 1)
        stop("Argument maf_max out of range")
    if (sum(rho_skatO<0 | rho_skatO>1) != 0)
        stop("Argument rho_skatO out of range")

    ### dimension check ###
    if(dim(aaf)[1] != 16569)
        stop("Argument aaf dimension out of range")
    if(dim(coverage)[1] != 16569)
        stop("Argument coverage dimension out of range")

    if(!setequal(dim(aaf),dim(coverage)))
        stop("The dimensions of aaf and coverage do not match")

    ###

    ID_samples <- colnames(aaf)
    coverage <- coverage[,ID_samples]
    pheno <- pheno[ID_samples,]
    region <- unique(region)
    if(!is.null(kins)){
        kins <- kins[ID_samples , ID_samples]
    }

    rownames(aaf) <- c(1:.mtLength)
    rownames(coverage) <- c(1:.mtLength)

    if(G_coding =="heter"){
        aaf_cat <- ((aaf >= thre_lower) & (aaf <= thre_upper))
        aaf_cat <- aaf_cat + 0

    }else if(G_coding=="heter2"){
        aaf_cat <- aaf
        aaf_cat[aaf < thre_lower]<-0
        aaf_cat[aaf > thre_upper]<-0

    }

    rownames(aaf_cat) <- rownames(aaf)
    colnames(aaf_cat) <- colnames(aaf)

    aaf_cat[ coverage < coverage.qc ] <- NA

    ## The reads of these loci are not reliable, so these loci should be removed
    loci_removed <- .mtLociNUMT   # c(301,302,310,316,3107,16182)
    loci_removed <- as.character(loci_removed)
    loci_left <- setdiff(rownames(aaf_cat), loci_removed)
    aaf_cat <- aaf_cat[loci_left,]

    # only include variants
    aaf_cat <- aaf_cat[apply(aaf_cat, 1, var, na.rm=T)>0,]

    if(heter_scale){
        aaf_cat <- t(apply(aaf_cat, 1, scale))
    }

    maf_heter<-rowSums(aaf_cat, na.rm = T)/(ncol(aaf_cat))
    maf_heter[maf_heter>0.5]<-1-maf_heter[maf_heter>0.5]
    aaf_cat <-aaf_cat[maf_heter<=maf_max,]


    # function to generate SNPInfo for seqMeta
    aaf_cat_SNPInfo<-as.data.frame(matrix(0, nrow(aaf_cat) , 2))
    colnames(aaf_cat_SNPInfo)<-c("Name" , "gene")
    aaf_cat_SNPInfo$Name<-rownames(aaf_cat)
    for(i in 1:nrow(aaf_cat_SNPInfo)){
        aaf_cat_SNPInfo$gene[i] <- .AA$Gene[.AA$Pos==as.numeric(aaf_cat_SNPInfo$Name[i])]
    }
    rownames(aaf_cat_SNPInfo) <- as.character(aaf_cat_SNPInfo$Name)

    # only include loci of specified region
    aaf_cat_region <- aaf_cat[rownames(aaf_cat)%in%as.character(region) , ]
    aaf_region <- aaf[rownames(aaf_cat_region),]
    aaf_cat_SNPInfo <- aaf_cat_SNPInfo[rownames(aaf_cat_region),]

    seqMeta_cov <- as.data.frame(pheno)

    # weights of the test
    wts1 <- 1

    aaf_cat_region_t<-t(aaf_cat_region)


    formula_null <- as.formula(paste0(trait, " ~ ", paste0(covars, collapse=" + ")))



        tryCatch({
            scores<-seqMeta::prepScores2(Z = aaf_cat_region_t, formula=formula_null, family = family,
                                SNPInfo = aaf_cat_SNPInfo, snpNames = "Name", kins=kins,
                                aggregateBy = "gene", sparse = TRUE, data = seqMeta_cov)

            burden_test1<-seqMeta::burdenMeta(scores, SNPInfo = aaf_cat_SNPInfo, wts = wts1,
                                     snpNames = "Name" ,aggregateBy = "gene",
                                     verbose = FALSE)

            skat_test1<-seqMeta::skatMeta(scores, SNPInfo = aaf_cat_SNPInfo, wts = wts1,
                                 snpNames = "Name" ,aggregateBy = "gene",
                                 verbose = FALSE)

            skatO_test1<-seqMeta::skatOMeta(scores, rho=rho_skatO,
                                   SNPInfo = aaf_cat_SNPInfo,
                                   burden.wts = wts1, skat.wts=wts1,
                                   snpNames = "Name" ,aggregateBy = "gene",
                                   verbose = FALSE, method="liu")

        }, error = function(e) {print(paste("error",sep=""))})

        beta_combined_burden1 <- burden_test1$beta
        se_combined_burden1 <- burden_test1$se
        p_combined_burden1 <- burden_test1$p

        Q_combined_skat1 <- skat_test1$Q
        p_combined_skat1 <- skat_test1$p

        p_combined_skatO1 <- skatO_test1$p

        # combine p values of tests of a gene/region by ACAT
        p_acat_all <- cbind(burden_test1$p,
                            skat_test1$p)
        Q_acat <- rowSums(tan((0.5-p_acat_all)*pi))/ncol(p_acat_all)
        p_combined_acat <- 0.5-(atan(Q_acat)/pi)

        ## output length, number of heteroplasmic loci in each gene,
        ## distribution of heteroplasmic loci w.r.t AAF range

            out_all_res <- cbind(beta_combined_burden1, se_combined_burden1, p_combined_burden1,
                                 Q_combined_skat1, p_combined_skat1,
                                 p_combined_skatO1,
                                 p_combined_acat)


        rownames(out_all_res) <- burden_test1$gene
        return(out_all_res)


}







