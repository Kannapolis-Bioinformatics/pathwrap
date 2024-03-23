#' Make TxDb object either from annotation file or
#' by loading from the annotation package
#'
#' Make_txdbobj uses the makeTxDbFromGFF to make
#' TxDb object from transcript annotations
#' available in gtf file in ref.dir or if ref.dir is NA
#' It makes txdb object fromt the annotaion package
#'
#' @param geneAnnotation : annotation file or if using package, path to sqlite
#' file in package
#' @param entity : the scientific name of the organism
#' @param outdir : directory to store output, Results is default
#' @import GenomicFeatures
#' @importFrom AnnotationDbi loadDb
#' @importFrom AnnotationDbi saveDb
#' @return txdb object that is returned
#'

make_txdbobj <-
    function(geneAnnotation,entity, outdir) {
        options(cache_size = NULL, synchronous = NULL)
        txdb <- try(loadDb(geneAnnotation), silent = TRUE)
        if (inherits(txdb, "TxDb")){
            message("Annotation generated from bioconductor package")
        }else{
            txdb <- makeTxDbFromGFF( 
                            file = geneAnnotation,
                            format = "auto",
                            chrominfo = NULL,
                            dataSource = NA,
                            organism = entity
            )
        }
        AnnotationDbi::saveDb(txdb,
            file = file.path(outdir, paste0(gsub(" ", "", entity), "_txdbobj"),
            fsep = .Platform$file.sep) )
        return(txdb)
}
