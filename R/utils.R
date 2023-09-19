#' function to make sure things are good
#'
#' @param ref.dir : directory for reference files
#' @param outdir : directory for result
#' @param entity : scientific name of species of interest
#' @param corenum : number of cores
#' @param compare : comparision to make for gage
#' @return list of different paths for result files
#' @importFrom stringr str_replace_all
#' @import Rsamtools
#' @importFrom methods is
#'

sanity_check <- function(ref.dir, outdir, entity, corenum, compare){
    # References
    # if only species name is given and both geneAnnotation and genome is NULL
    if (is.na(ref.dir)) {
    data(anntpkglist, package = "pathviewwrap", envir = environment())
    ref_info <- anntpkglist
    species_no <- which(ref_info$species == entity)
    annotate_pkg <- ref_info$annotation[species_no]
    genome_pkg <- ref_info$genome[species_no]
        # annotation pkg installation
        pkg.on <- requireNamespace(annotate_pkg,
            lib.loc = .libPaths()[1],
            quietly = TRUE
        )
        if (!pkg.on) {
            message(paste0("
        Intall the required package with the following command,
        > BiocManager::install('", annotate_pkg, "'
                        ,force = TRUE,
                        lib.loc = .libPaths()[1]  )"
            , collapse = ""))
            return(invisible(NULL)) # needs reevaluation
        }
        geneAnnotation <- file.path(
            .libPaths()[1], annotate_pkg,
            "extdata", paste0(annotate_pkg, ".sqlite")
        )
        # genome file installation
        genomeFile <- genome_pkg
        pkg.on <- requireNamespace(genome_pkg,
            lib.loc = .libPaths()[1],
            quietly = TRUE
        )
        if (!pkg.on) {
            message(paste0(
                "
        Intall the required package with the following command,
        >  BiocManager::install('",
                genome_pkg, "',force = TRUE,
                        lib.loc = .libPaths()[1] )",collapse = ""
            ))
            return(invisible(NULL))
        }
    } else {
        geneAnnotation <- list.files(ref.dir, ".gtf$|.gff$", 
                                     full.names = TRUE)[1]
        genomeFile <- list.files(ref.dir, ".fa$|.fna$|.fa.gz",
            full.names = TRUE)[1]
       
        if (summary(file(genomeFile))$class == "gzfile") {
           a <- try({
                system(paste0("gunzip -k ", genomeFile))
                genomeFile <- str_remove(pattern = ".gz$", genomeFile)
            }, silent = TRUE)
           loadError <- (is(a, "try-error") | is(a, "error"))
           if (loadError == TRUE) {
               message("GenomeFile should be bgzip file not gzipped")
               return(invisible(x = NULL ))
           }
        }
    }
    
    return(c(genomeFile, geneAnnotation))
}




#' Function to clean the reference directories/packages
#' @param ref.dir : directory for reference files
onexistcleanup <- function(ref.dir, entity){
    if (is.na(ref.dir)) {
        data(anntpkglist, package = "pathviewwrap", envir = environment())
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
                                  annotate_pkg, "extdata", sqlite.md5
        ))) {
            unlink(file.path(.libPaths()[1],
                             annotate_pkg, "extdata", sqlite.md5
            ))}
        if (file.exists(file.path(.libPaths()[1],
                                  annotate_pkg, "extdata",
                                  sqlite.SpliceSites.txt.md5
        ))) {
            unlink(file.path( .libPaths()[1], annotate_pkg,
                              "extdata", sqlite.SpliceSites.txt.md5
            ))}
        if (file.exists(file.path(.libPaths()[1], annotate_pkg,
                                  "extdata", sqlite.SpliceSites.txt
        ))) {
            unlink(file.path(
                .libPaths()[1], annotate_pkg,
                "extdata", sqlite.SpliceSites.txt
            ))
            }
    }
    }
