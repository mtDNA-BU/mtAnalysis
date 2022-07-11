# function to run association analysis of heteroplasmies





#' mtAnalysis function
#'
#' @param aaf test1
#' @param family
#' @param Methods
#' @param pheno
#' @param trait
#' @param covars test2
#' @param G_coding
#' @param region
#' @param aaf_cat_SNPInfo
#' @param kins
#'
#' @return
#' @export
#'
#' @examples
mtAnalysis<-function(aaf, family, Methods, pheno, trait, covars, G_coding,
                     region, aaf_cat_SNPInfo=aaf_cat_SNPInfo,
                     kins=NULL){

    if(G_coding =="heter"){
        aaf_cat <- ((aaf >= thre_spec_lower) & (aaf <= thre_spec_upper))
        aaf_cat <- aaf_cat + 0

    }else if(G_coding=="heter2"){
        aaf_cat <- aaf
        aaf_cat[aaf < thre_spec_lower]<-0
        aaf_cat[aaf > thre_spec_upper]<-0

    }

    rownames(aaf_cat) <- rownames(aaf)
    colnames(aaf_cat) <- colnames(aaf)

    # only include variants
    aaf_cat <- aaf_cat[apply(aaf_cat, 1, var, na.rm=T)>0,]

    # function to generate SNPInfo for seqMeta
    aaf_cat_SNPInfo<-as.data.frame(matrix(0, nrow(aaf_cat) , 2))
    colnames(aaf_cat_SNPInfo)<-c("Name" , "gene")
    aaf_cat_SNPInfo$Name<-rownames(aaf_cat)
    for(i in 1:nrow(aaf_cat_SNPInfo)){
        aaf_cat_SNPInfo$gene[i] <- AA$Gene[AA$Pos==as.numeric(aaf_cat_SNPInfo$Name[i])]
    }
    rownames(aaf_cat_SNPInfo) <- as.character(aaf_cat_SNPInfo$Name)

    # only include loci of specified region
    aaf_cat_region <- aaf_cat[rownames(aaf_cat)%in%as.character(region) , ]
    aaf_region <- aaf[rownames(aaf_cat_region),]
    aaf_cat_SNPInfo <- aaf_cat_SNPInfo[rownames(aaf_cat_region),]
    HMF_bar <- HMF_bar[rownames(aaf_cat_region)]
    mt_test_names2 <- unique(aaf_cat_SNPInfo$gene)

    seqMeta_cov <- as.data.frame(pheno)
    num_genes <- length(unique(aaf_cat_SNPInfo$gene))

    aaf_cat_region_count <- ((aaf_region>=thre_count_lower) & (aaf_region<=thre_count_upper))
    aaf_cat_region_count <- aaf_cat_region_count + 0

    # weights of the test
    wts1 <- 1

    aaf_cat_region_t<-t(aaf_cat_region)


    formula_null <- as.formula(paste0("y ~ ", covars))

    if(Methods=="common"){

        seqMeta_cov$y <- pheno[,trait]
        tryCatch({
            scores<-seqMeta::prepScores2(Z = aaf_cat_region_t, formula=formula_null, family = family,
                                SNPInfo = aaf_cat_SNPInfo, snpNames = "Name",
                                aggregateBy = "gene", sparse = TRUE, data = seqMeta_cov)

            burden_test1<-seqMeta::burdenMeta(scores, SNPInfo = aaf_cat_SNPInfo, wts = wts1,
                                     snpNames = "Name" ,aggregateBy = "gene",
                                     verbose = FALSE)

            skat_test1<-seqMeta::skatMeta(scores, SNPInfo = aaf_cat_SNPInfo, wts = wts1,
                                 snpNames = "Name" ,aggregateBy = "gene",
                                 verbose = FALSE)

            skatO_test1<-seqMeta::skatOMeta(scores, rho=c(0,1),
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
        if(G_coding=="heter"){
            len_genes <- rep(NA, num_genes)
            for(i in seq_len(num_genes)){
                len_genes[i] <- sum(AA$gene==mt_test_names2[i])
            }
            n_mut_genes <- table(aaf_cat_SNPInfo$gene)
            n_mut_genes <- n_mut_genes[mt_test_names2]
            n_mut_genes <- as.vector(n_mut_genes)

            aaf_cat_dot_pro <- aaf_region*aaf_cat_region_count
            aaf_mut_upper <- apply(aaf_cat_dot_pro, 1, FUN=max, na.rm=T)
            aaf_mut_lower <- rep(NA, nrow(aaf_cat_dot_pro))
            for(i in 1:nrow(aaf_cat_dot_pro)){
                aaf_cat_dot_pro_row <- aaf_cat_dot_pro[i,]
                aaf_cat_dot_pro_row <- aaf_cat_dot_pro_row[aaf_cat_dot_pro_row>0]
                aaf_mut_lower[i] <- min(aaf_cat_dot_pro_row, na.rm=T)
            }
            n_mut_0325 <- rep(NA, num_genes)
            n_mut_0350 <- rep(NA, num_genes)
            n_mut_0375 <- rep(NA, num_genes)
            n_mut_0397 <- rep(NA, num_genes)
            n_mut_5097 <- rep(NA, num_genes)
            for(i in seq_len(num_genes)){
                aaf_mut_upper_gene <- aaf_mut_upper[aaf_cat_SNPInfo$gene==mt_test_names2[i]]
                aaf_mut_lower_gene <- aaf_mut_lower[aaf_cat_SNPInfo$gene==mt_test_names2[i]]
                n_mut_0325[i] <- sum((aaf_mut_upper_gene<=0.25), na.rm=T)
                n_mut_0350[i] <- sum((aaf_mut_upper_gene>0.25) & (aaf_mut_upper_gene<=0.5), na.rm=T)
                n_mut_0375[i] <- sum((aaf_mut_upper_gene>0.5) & (aaf_mut_upper_gene<=0.75), na.rm=T)
                n_mut_0397[i] <- sum((aaf_mut_lower_gene<=0.5) & (aaf_mut_upper_gene>0.75), na.rm=T)
                n_mut_5097[i] <- sum((aaf_mut_lower_gene>0.5) & (aaf_mut_upper_gene>0.75), na.rm=T)
            }

            out_all_res <- cbind(beta_combined_burden1, se_combined_burden1, p_combined_burden1,
                                 Q_combined_skat1, p_combined_skat1,
                                 p_combined_skatO1,
                                 p_combined_acat,
                                 len_genes, n_mut_genes, n_mut_0325, n_mut_0350,
                                 n_mut_0375, n_mut_0397, n_mut_5097)
        }else{
            out_all_res <- cbind(beta_combined_burden1, se_combined_burden1, p_combined_burden1,
                                 Q_combined_skat1, p_combined_skat1,
                                 p_combined_skatO1,
                                 p_combined_acat)
        }

        rownames(out_all_res) <- burden_test1$gene
        return(out_all_res)

    }
}
























