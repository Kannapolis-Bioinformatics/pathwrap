#' Make TxDb object either from annotation file or
#' by loading from the annotation package
#'
#' Make_txdbobj uses the makeTxDbFromGFF to make
#' TxDb object from transcript annotations
#' available in gtf file in ref.dir or if ref.dir is NA
#' It makes txdb object fromt the annotaion package
#'
#' @param geneAnnotation : annotation file or package
#' @param corenum : the number of cores
#' @param genomeFile : genomeFile or package
#' @param entity  : the scientific name of the organism
#' @param outdir : directory to store output, Results is default
#'
#' @import GenomicFeatures
#' @import Rsamtools
#' @importFrom AnnotationDbi loadDb
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb renameSeqlevels
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics width
#' @return txdb object that is returned
#'

make_txdbobj <- function(geneAnnotation, corenum, genomeFile, entity, outdir) {
    options(cache_size = NULL, synchronous = NULL)
    txdb <- try(loadDb(geneAnnotation), silent = TRUE)
    cl2 <- makeCluster(corenum)
    # if (class(txdb)==  "TxDb"){
    if (is(txdb, "TxDb")) {
        if (!grepl("chr", seqlevels(txdb)[1])) { # check if this is necessary
            newSeqNames <- paste("Chr", seqlevels(txdb), sep = "")
            names(newSeqNames) <- seqlevels(txdb)
            txdb <- renameSeqlevels(txdb, newSeqNames)
            # seqlevels(txdb)
        }
        closeAllConnections()
    } else {
        chrLen <- Rsamtools::scanFaIndex(genomeFile)
        chrominfo <- data.frame(
            chrom = as.character(seqnames(chrLen)),
            length = width(chrLen),
            is_circular = rep(FALSE, length(chrLen))
        )
        txdb <- makeTxDbFromGFF(
            file = geneAnnotation, format = "gtf",
            chrominfo = chrominfo,
            dataSource = "Ensembl",
            organism = entity
        )
    }
    AnnotationDbi::saveDb(txdb,
        file = paste0(outdir, "/", gsub(" ", "", entity), "_txdbobj")
    )
    return(txdb)
}
