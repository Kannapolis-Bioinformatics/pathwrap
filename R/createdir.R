#' create the directories for saving the result
#'
#' @param outdir where the directory and its subfolders will be created
#' @param startover if TRUE, previous results in outdir (if any) are cleared
#' @return the resolved outdir path (invisible NULL if the user declines)
createdir <- function(outdir, startover) {
    
    #if (is.null(outdir)) {
    #    outdir <- file.path(getwd(), "results")
    #}
    if (!file.exists(outdir)) {
        ok <- dir.create(outdir, recursive = TRUE)
        if (!ok || !file.exists(outdir)) {
            stop("Failed to create output directory: ", outdir,
                 " - check that the parent path is valid and writable.")
        }
    }
    message("This is outdir: ", outdir)
    
    outdir_has_content <- file.exists(outdir) &&
        length(list.files(outdir, all.files = TRUE, no.. = TRUE)) > 0
    
    if (startover) {
        if (outdir_has_content) {
            if (!interactive() || getOption("pathwrap.noninteractive", FALSE)) {
                # Docker / Rscript / CI - proceed silently
                message("Non-interactive session: clearing existing contents of ", outdir)
                unlink(outdir, recursive = TRUE)
            } else {
                answer <- readline("outdir is not empty. Delete everything in it? (yes/no): ")
                if (tolower(answer) != "yes") {
                    stop("startover=TRUE was declined. Set startover=FALSE and keep the ",
                         "existing results, or confirm deletion to start clean.")
                }
                unlink(outdir, recursive = TRUE)
            }
            ok <- dir.create(outdir, recursive = TRUE)
            if (!ok || !file.exists(outdir)) {
                stop("Failed to recreate output directory: ", outdir)
            }
        }
        # outdir was already empty - nothing to delete, nothing to confirm
    }
    
    if (!startover) {
        if (outdir_has_content) {
            if (!interactive() || getOption("pathwrap.noninteractive", FALSE)) {
                message("Existing results found in ", outdir, ".")
                message("Non-interactive session: proceeding without manual confirmation.")
                message("Parameters will still be checked automatically against the previous run's manifest.")
            } else {
                ans <- readline("Existing results found. Is this the same run with the same parameters as before? (yes/no): ")
                if (substr(tolower(ans), 1, 1) == "n") {
                    stop("Make sure the parameters are the same as the previous run, ",
                         "or answer yes above to resume.")
                }
            }
        }
    }
    
    message(paste0("The results will be organized in ", outdir, collapse = "|"))
    
    checkcretdir <- function(parentname, dirname) {
        full <- file.path(parentname, dirname, fsep = .Platform$file.sep)
        if (!file.exists(full)) dir.create(full, recursive = TRUE)
        invisible(full)
    }
    
    folder_to_create <- list("fastqc_results", "fastp_results",
                             "gage_results", "differential_analysis", "aligned_bam")
    trim_dir <- list("fastp_failed", "unpaired")
    diff_dir <- list("DESeq2", "edgeR")
    pathway_types <- list("KEGG", "GO", "kegg_csets", "combined_analysis_kegg")
    kegg_types <- list("signalling", "metabolism", "disease", "signmet")
    go_types <- list("biological_process", "molecular_function", "cellular_component")
    
    lapply(folder_to_create, checkcretdir, parentname = outdir)
    lapply(trim_dir, checkcretdir, parentname = file.path(outdir, "fastp_results", fsep = .Platform$file.sep))
    lapply(diff_dir, checkcretdir, parentname = file.path(outdir, "differential_analysis", fsep = .Platform$file.sep))
    lapply(pathway_types, checkcretdir, parentname = file.path(outdir, "gage_results", fsep = .Platform$file.sep))
    lapply(kegg_types, checkcretdir, parentname = file.path(outdir, "gage_results", "KEGG", fsep = .Platform$file.sep))
    lapply(go_types, checkcretdir, parentname = file.path(outdir, "gage_results", "GO", fsep = .Platform$file.sep))
    #lapply(list("user_pathways"), checkcretdir, parentname = file.path(outdir, "gage_results", fsep = .Platform$file.sep))
    return(outdir)
}