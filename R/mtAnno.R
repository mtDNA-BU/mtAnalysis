#' mtAnno function
#'
#' Annotate Given Alternative Alleles of mtDNA Loci
#'
#' @param anno a data frame provided by the user. It has two columns:
#' mtDNA loci positions ("pos") and
#' alternative alleles ("alleles") to be annotated.
#' @param annot.select types of annotation scores to be generated.
#' @param path the path to the directory of the output of annotation file.
#' If not provided, the annotation
#' file will output to the current working directory.
#' @param study A string of study names. Default is "Study".
#' @return A .csv file containing annotated alternative alleles with
#' corresponding mtDNA loci positions.
#' @export
#' @examples
#'
#'\dontrun{
#' ## Read input data
#' allele_file <- "allele.csv"
#' freq_file   <- "freq.csv"
#'
#' allele <- as.matrix( read.csv(file = allele_file, sep = ",", header=FALSE) )
#' freq   <- as.matrix( read.csv(file = freq_file, sep = ",", header=FALSE) )
#'
#' aaf =  mtAAF ( allele, freq)
#'
#' ## Annotate Given Alternative Alleles of mtDNA Loci
#' mtAnno(aaf)
#'
#' ##
#' mtAnno(aaf, annot.select=c("ref","Gene") )
#'}

mtAnno <- function(anno,
                   annot.select=c("Pos", "ref", "Gene", "TypeMutation",
                                  "MissensMutation", "CodonPosition",
                                  "ProteinDomain", "dbSNP_150_id", "PolyPhen2",
                                  "PolyPhen2_score", "SIFT","SIFT_score",
                                  "CADD", "CADD_score", "CADD_phred_score"),
                   path="./", study="Study"){

    anno$alleles <- as.character(anno$alleles)
    anno_score <- NULL

    ## annotate the alleles based on the annotation files
    for(i in  1:nrow(anno)){
        if(anno$alleles[i] == "A"){
            score <- .AA[.AA$Pos == anno$pos[i], annot.select]
        }
        else if(anno$alleles[i] == "T"){
            score <- .TT[.TT$Pos == anno$pos[i], annot.select]
        }
        else if(anno$alleles[i] == "G"){
            score <- .GG[.GG$Pos == anno$pos[i], annot.select]
        }
        else if(anno$alleles[i] == "C"){
            score <- .CC[.CC$Pos == anno$pos[i], annot.select]
        }
        anno_score <- rbind(anno_score, score)
    }

    output <- cbind(anno, anno_score)
    write.csv(output, file=paste0(path, study, "_mtAnno.csv"), row.names = F)

}





