#' function to make sure things are good
#'
#' @param ref.dir : directory for reference files
#' @param outdir : directory for result
#' @param entity : scientific name of species of interest
#' @param corenum : number of cores
#' @param compare : comparision to make for gage
#' @param pos : environment specifier
#' @param rerun if FALSE the previously complete step will not be rerunned, if
#' TRUE analysis starts from first step
#' @return list of different paths for result files
#' @importFrom stringr str_replace_all
#' @import BiocManager
#' @import Rsamtools
#' @importFrom methods is
#'

sanity_check <- function(ref.dir, outdir, pos = 1, entity, corenum, compare,
                        rerun) {
    if (file.exists(outdir) & rerun == FALSE) {
        unlink(outdir, recursive = TRUE)
    }

    if (!file.exists(outdir)) {
        # default output file
        dir.create(outdir)
    }
    result.dir <- outdir
    message(paste0("The results will be organized in ", result.dir))
    # setwd(outdir)

    # check and create dir for organizing results
    checkcretdir <- function(parentname, dirname) {
        if (!file.exists(file.path(parentname, dirname))) {
            dir.create(file.path(parentname, dirname))
        }
        assign(dirname,
            value = file.path(parentname, dirname),
            envir = as.environment(pos)
        ) # .GlobalEnv)#environment())
    }

    folder_to_create <- list(
        "fastqc_results", "fastp_results", "gage_results",
        "differential_analysis", "aligned_bam"
    )
    trim_dir <- list("fastp_log", "unpaired")
    diff_dir <- list("DESeq2", "edgeR")
    pathway_types <- list("KEGG", "GO")
    kegg_types <- list("signalling", "metabolism", "disease", "sig_n_met")
    go_types <- list(
        "biological_process", "molecular_function",
        "cellular_component"
    )
    lapply(folder_to_create, checkcretdir, parentname = result.dir)
    lapply(trim_dir, checkcretdir, parentname = file.path(
        result.dir,
        "fastp_results"
    ))
    lapply(diff_dir, checkcretdir, parentname = file.path(
        result.dir,
        "differential_analysis"
    ))
    lapply(pathway_types, checkcretdir, parentname = file.path(
        result.dir,
        "gage_results"
    ))
    lapply(kegg_types, checkcretdir, parentname = file.path(
        result.dir,
        "gage_results", "KEGG"
    ))
    lapply(go_types, checkcretdir, parentname = file.path(
        result.dir,
        "gage_results", "GO"
    ))

    # just to make sure rest of codes are same
    qc.dir <- fastqc_results
    diff.dir <- differential_analysis
    trim.dir <- fastp_results
    gage.dir <- gage_results
    trim.log <- fastp_log
    edger.dir <- edgeR
    message(edger.dir)
    deseq2.dir <- DESeq2
    kegg.dir <- KEGG
    go.dir <- GO

    # References
    # if only species name is given and both geneAnnotation and genome is NULL
    if (is.na(ref.dir)) {
        data(anntpkglist, package = "pathviewwrap", envir = environment())
        ref_info <- anntpkglist

        species_no <- which(ref_info$species == entity)
        annotate_pkg <- ref_info$annotation[species_no]
        genome_pkg <- ref_info$genome[species_no]

        # (set of genome and annotation pkg come from developers list)
        #
        sqlite.md5 <- paste0(annotate_pkg, ".sqlite.md5")
        sqlite.SpliceSites.txt.md5 <- paste0(
            annotate_pkg,
            ".sqlite.SpliceSites.txt.md5"
        )
        sqlite.SpliceSites.txt <- paste0(
            annotate_pkg,
            ".sqlite.SpliceSites.txt"
        )
        if (file.exists(file.path(
            .libPaths()[1],
            annotate_pkg, "extdata", sqlite.md5
        ))) {
            unlink(file.path(
                .libPaths()[1],
                annotate_pkg, "extdata", sqlite.md5
            ))
        }
        if (file.exists(file.path(
            .libPaths()[1],
            annotate_pkg, "extdata",
            sqlite.SpliceSites.txt.md5
        ))) {
            unlink(file.path(
                .libPaths()[1], annotate_pkg,
                "extdata", sqlite.SpliceSites.txt.md5
            ))
        }
        if (file.exists(file.path(
            .libPaths()[1], annotate_pkg,
            "extdata", sqlite.SpliceSites.txt
        ))) {
            unlink(file.path(
                .libPaths()[1], annotate_pkg,
                "extdata", sqlite.SpliceSites.txt
            ))
        }

        # annotation pkg installation
        pkg.on <- requireNamespace(annotate_pkg,
            lib.loc = .libPaths()[1],
            quietly = TRUE)
        if (!pkg.on) {
            message(paste0("
        Intall the required package with the following command,
        > BiocManager::install('", annotate_pkg, "'
                        ,force = TRUE, 
                        lib.loc = .libPaths()[1]  )"))
            return(invisible(NULL)) # needs reevaluation
            # pkg.on = requireNamespace(annotate_pkg,  lib.loc = .libPaths()[1])
            # if (!pkg.on)
            #  stop(paste("Fail to install/load gene annotation package
            # ",annotate_pkg, "!", sep = ""))
        }
        geneAnnotation <- file.path(
            .libPaths()[1], annotate_pkg,
            "extdata", paste0(annotate_pkg, ".sqlite")
        )

        # genome file installation
        genomeFile <- genome_pkg
        pkg.on <- requireNamespace(genome_pkg,
            lib.loc = .libPaths()[1],
            quietly = TRUE)
        if (!pkg.on) {
            message(paste0("
        Intall the required package with the following command,
        >  BiocManager::install('",
                genome_pkg, "',force = TRUE,
                        lib.loc = .libPaths()[1] )"
            ))
          return(invisible(NULL))
        }
    } else {
        genomeFile <- list.files(ref.dir, ".fa$|.fna$|.fa.gz",
            full.names = TRUE
        )[1]
        # unzipping .gz file because both scanFaIndex and qAlign
        # do not work with gzip' ed file, require bgzip file
        if (summary(file(genomeFile))$class == "gzfile") {
            system(paste0("gunzip -k ", genomeFile))
            genomeFile <- str_remove(pattern = ".gz$", genomeFile)
        }
        geneAnnotation <- list.files(ref.dir, ".gtf$|.gff$", full.names = TRUE)
        ## could be changed to include gtf, gff etc,check with quasR package
        message(geneAnnotation)
    }
    message(paste0(
        "this is directory lists,
        qc.dir :", qc.dir, "\n", 
        "trim.dir :", trim.dir, "\n", 
        "genomeFile :", genomeFile, "\n", 
        "geneAnnotation :", geneAnnotation, "\n", 
        "deseq2.dir :", deseq2.dir, "\n", 
        "edger.dir :", edger.dir, "\n", 
        "gage.dir : ", gage.dir, "\n"
    ))

    return(c(
        qc.dir, trim.dir, genomeFile, geneAnnotation, deseq2.dir,
        edger.dir, gage.dir
    ))
}
