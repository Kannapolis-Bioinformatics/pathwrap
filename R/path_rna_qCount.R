#' Count the number of genes expression and store count in RDS file
#'
#' generate the count of gene and if gene id is ensembl, it converts it to
#' entrez to match with genes names for gene set analysis.
#' the generated table with genes in row and counts in columns are stored as a
#' RDS file that can be loaded into R for future analysis.
#'
#' @param aligned_proj_lst :the name of R object saved or generated after all alignments
#' @param corenum : the number of cores used for alignment
#' @param outdir : the directory in which result is stored
#' @param references : the references which is used for alignment and counting
#' @param entity : the scientific name of the species of interest
#'
#' @import stringr
#' @importFrom QuasR qCount
#' @importFrom Rsamtools scanFaIndex
#' @import GenomicFeatures
#' @import gage
#' @import pathview
#' @import parallel
#'
#' @return count of genes
#'

run_qCount <- function(aligned_proj_lst, corenum, outdir ,entity, references) {
    ##
    # for mapping
    # library(Rsamtools) #scanFaIndex
    if(!file.exists(file.path(outdir,"aligned_bam", "combinedcount.trimmed.RDS",
                             fsep = .Platform$file.sep))){
      txdb <-make_txdbobj(references$geneAnnotation,entity, outdir)
      
      cnts_all <- c()
      for (aligned_proj in aligned_proj_lst){
        cl2 <- makeCluster(corenum)
        geneLevels <- QuasR::qCount(aligned_proj, txdb,
          reportLevel = "gene",collapseBySample=FALSE,
          clObj = cl2
        )
        stopCluster(cl2)
      #####################
      # post processing for count
        oldcolnames <-colnames(geneLevels)
        cnts <- as.matrix(geneLevels[, -1])
        colnames(cnts)<- oldcolnames[-1]
        
        cnts_all <- cbind(cnts_all,cnts)
      }
      #cnts_all <- do.call(cbind, cnts_all)
      
      
      saveRDS(cnts_all, file.path(outdir,"aligned_bam", "combinedcount.trimmed.RDS",
                              fsep = .Platform$file.sep))
    } else {
      cnts_all<-readRDS(file.path(outdir,"aligned_bam", "combinedcount.trimmed.RDS",
                          fsep = .Platform$file.sep))
    }
    message("the number of genes counted is ", dim(cnts_all)[1])
    return(cnts_all)
}
