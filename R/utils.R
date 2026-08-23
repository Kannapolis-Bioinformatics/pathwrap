#' function to make sure things are good
#'
#' @param ref.dir : directory for reference files
#' @param entity : scientific name of species of interest
#' @return list of different paths for result files
#' @importFrom stringr str_replace_all
#' @import Rsamtools
#' @importFrom methods is
#'

check_references <- function(ref.dir, entity){
    # References
    # if only species name is given and both geneAnnotation and genome is NULL
    if (is.null(ref.dir)) {
        data(anntpkglist, package = "pathwrap", envir = environment())
        ref_info <- anntpkglist
        species_no <- which(ref_info$species == entity)
        annotate_pkg <- ref_info$annotation[species_no]
        genome_pkg <- ref_info$genome[species_no]
        # annotation pkg installation
        pkg.on <- requireNamespace(annotate_pkg,lib.loc = .libPaths()[1],
                                   quietly = TRUE)
        if (!pkg.on) {
            try.on  <- try( BiocManager::install( annotate_pkg, 
                                                  force = TRUE,
                                                  lib = .libPaths()[1] , 
                                                  ask= FALSE )
            )
            if (inherits(try.on, "try-error")){
                
                message(paste0("
                Intall the required package with the following command,
                > BiocManager::install('", annotate_pkg, "'
                        ,force = TRUE,
                        lib= .libPaths()[1]  )", collapse = ""))
                return(invisible(NULL)) 
            }
        }
        geneAnnotation <- file.path(
            .libPaths()[1], annotate_pkg, "extdata", 
            paste0(annotate_pkg, ".sqlite"), fsep= .Platform$file.sep)
        genomeFile <- genome_pkg # genome file installation
        pkg.on <- requireNamespace(genome_pkg,
                                   lib.loc = .libPaths()[1], quietly = TRUE)
        if (!pkg.on) {
            message(paste0(
                "Intall the required package with the following command,
        >  BiocManager::install('",
                genome_pkg, "',force = TRUE,
                        lib = .libPaths()[1] )",collapse = ""))
            return(invisible(NULL)) }
    } else {
            geneAnnotation <- list.files(ref.dir, "\\.gtf$|\\.gff$", 
                                         full.names = TRUE)[1]
            genomeFile <- list.files(ref.dir, "\\.fa$|\\.fna$|\\.fa.gz$",
                                     full.names = TRUE)[1]
            
                con <- file(genomeFile)
                on.exit(close(con), add = TRUE)
                if (summary(con)$class == "gzfile"){
                    a <- try({
                        exit_code<-system(paste0("gunzip -k ", genomeFile))
                        if (exit_code != 0) stop("gunzip exited with code ", 
                                                 exit_code)
                        genomeFile <- str_remove(pattern = ".gz$", genomeFile)
                    }, silent = TRUE)
                    loadError <- (is(a, "try-error") | is(a, "error"))
                    if (loadError == TRUE) {
                        message("GenomeFile should be bgzip file not gzipped")
                        return(invisible(x = NULL ))
                    }
                }
                
        }
    return(list("genomeFile" = genomeFile,"geneAnnotation"= geneAnnotation))
}

#' Function to clean the reference directories/packages


#' Function to clean the reference directories/packages
#' @param entity : scientific name of the organism
#' @param ref.dir : directory for reference filesi
#' @return return message about clean up completion
#' @export
onexistcleanup <- function(ref.dir, entity){
    if (is.null(ref.dir)) {
        data(anntpkglist, package = "pathwrap", envir = environment())
        ref_info <- anntpkglist
        species_no <- which(ref_info$species == entity)
        annotate_pkg <- ref_info$annotation[species_no]
        genome_pkg <- ref_info$genome[species_no]
        # (set of genome and annotation pkg come from developers list)
        sqlite.md5 <- paste0(annotate_pkg, ".sqlite.md5")
        sqlite.SpliceSites.txt.md5 <- paste0(annotate_pkg,
                                             ".sqlite.SpliceSites.txt.md5")
        sqlite.SpliceSites.txt <- paste0(annotate_pkg,
                                         ".sqlite.SpliceSites.txt")
        if (file.exists(file.path(.libPaths()[1],
                                  annotate_pkg, "extdata", sqlite.md5, fsep = .Platform$file.sep
        ))) {
            unlink(file.path(.libPaths()[1],
                             annotate_pkg, "extdata", sqlite.md5,fsep = .Platform$file.sep
            ))}
        if (file.exists(file.path(.libPaths()[1],
                                  annotate_pkg, "extdata",
                                  sqlite.SpliceSites.txt.md5, fsep = .Platform$file.sep
        ))) {
            unlink(file.path( .libPaths()[1], annotate_pkg,
                              "extdata", sqlite.SpliceSites.txt.md5, fsep = .Platform$file.sep
            ))}
        if (file.exists(file.path(.libPaths()[1], annotate_pkg,
                                  "extdata", sqlite.SpliceSites.txt, fsep = .Platform$file.sep
        ))) {
            unlink(file.path(
                .libPaths()[1], annotate_pkg,
                "extdata", sqlite.SpliceSites.txt, fsep = .Platform$file.sep
            ))
        }
    }
    return ("package clean up complete.")
}