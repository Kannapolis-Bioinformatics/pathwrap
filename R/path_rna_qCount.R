#' Count the number of genes expression and store count in RDS file
#'
#' generate the count of gene and if gene id is ensembl, it converts it to
#' entrez to match with genes names for gene set analysis.
#' the generated table with genes in row and counts in columns are stored as a
#' RDS file that can be loaded into R for future analysis.
#'
#' @param aligned_proj :the name of R object saved or generated during alignment
#' @param corenum : the number of cores used for alignment
#' @param result.dir : the directory in which result is stored
#' @param txdb : the txdb object which is used for counting
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

run_qCount <- function(aligned_proj, corenum, result.dir, txdb, entity) {
    ##
    # for mapping
    # library(Rsamtools) #scanFaIndex
    cl2 <- makeCluster(corenum)
    geneLevels <- QuasR::qCount(aligned_proj, txdb,
        reportLevel = "gene",
        clObj = cl2
    )

    #####################
    # post processing for count
    cnts <- geneLevels[, -1]
    kegg.gs.species <- kegg.gsets(entity)
    orgcode <- kegg.species.code(entity)
    data(bods, package = "gage", envir = environment())
    # if(!all(rownames(cnts)%in% unlist(unname(kegg.gs.species$kg.sets))))
    # { #check if the use of "all" is appropriate
    if (sum(rownames(cnts) %in% unlist(unname(kegg.gs.species$kg.sets))) < 10) {
        rownames(cnts) <- str_remove(rownames(cnts), "\\.[0-9]+$")
        cnts <- mol.sum(cnts,
            id.map = "ENSEMBL",
            gene.annotpkg = bods[which(bods[, 3] == orgcode)]
        )
        # converting to entrez # what if gene id is not ensembl and what if
        # arabidopsis thaliana id.map might be ath or else thing
    }
    saveRDS(cnts, file.path(result.dir, "combinedcount.trimmed.RDS"))
    return(cnts)
}
