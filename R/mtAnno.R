#' mtAnno function
#'
#' Annotate mtDNA alternative alleles
#'
#' @param anno data frame provided by the user. It contains alleles need to be annotated.
#' It has two columns: loci positions ("pos") and "alleles" to be annotated.
#' @param annot.select types of annotation scores based on user's choice
#' @param path the path of the output annotation file. If not provided, the annotation file will
#' output to the current working directory
#' @param study A string of study names. Default is Study.
#' @return Annotation of alleles stored as a .csv file to the path provided by users
#' @export
#' @examples

mtAnno <- function(anno , annot.select=c("Pos","ref","Gene","TypeMutation","MissensMutation",
                                        "CodonPosition","ProteinDomain","dbSNP_150_id","PolyPhen2",
                                        "PolyPhen2_score","SIFT","SIFT_score", "CADD","CADD_score",
                                        "CADD_phred_score"), 
                   path="./",study="Study"){
  anno$alleles<-as.character(anno$alleles)
  anno_score<-NULL
  # annotate the alleles based on the annotation files
  for(i in 1:nrow(anno)){
    if(anno$alleles[i]=="A"){
      score<-.AA[.AA$Pos==anno$pos[i] , annot.select]
    }
    else if(anno$alleles[i]=="T"){
      score<-.TT[.TT$Pos==anno$pos[i] , annot.select]
    }
    else if(anno$alleles[i]=="G"){
      score<-.GG[.GG$Pos==anno$pos[i] , annot.select]
    }
    else if(anno$alleles[i]=="C"){
      score<-.CC[.CC$Pos==anno$pos[i] , annot.select]
    }
    anno_score<-rbind(anno_score , score)
  }
  
  # combine the input file with the annotation results
  output<-cbind(anno , anno_score)
  
  # output the annotated data
  write.csv(output , file=paste0(path , study ,"_mtAnno.csv"))
  
}





