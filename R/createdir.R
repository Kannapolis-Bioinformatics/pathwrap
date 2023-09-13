#' create the directories for saving the result
#'
#' @param entity : scientific name of the organism
#' @param rerun : if the run is again?
#' @param keep_tmp : keep bam files resulting from alignment or not
#' @param pos enviromental 
#' @param outdir where the directory and its subfolder will be created
#' @return list of folder names

createdir <- function(pos =1, outdir, entity, rerun, keep_tmp) {
    on.exit(closeAllConnections())
    aligned_bam <- NA
    if (keep_tmp == FALSE) {
        message("deleting aligned bam files, bam file index and log files")
        unlink(list.files(file.path(outdir, "aligned_bam"),
                        pattern = ".bam$|.bai$", full.names = TRUE
        ))
    }
    if (file.exists(outdir) & rerun == FALSE) {
        unlink(outdir, recursive = TRUE)
    }
    if (!file.exists(outdir)) {
        # default output file
        dir.create(outdir)
    }
    result.dir <- outdir
    message(paste0("The results will be organized in ", result.dir, 
        collapse = "|"))
    
    # check and create dir for organizing results
    checkcretdir <- function(parentname, dirname) {
        if (!file.exists(file.path(parentname, dirname))) {
            dir.create(file.path(parentname, dirname))
        }
        assign(dirname,
            value = file.path(parentname, dirname),
            envir = as.environment(pos)
        )} # .GlobalEnv)#environment())}
    folder_to_create <- list(
        "fastqc_results", "fastp_results", "gage_results",
        "differential_analysis", "aligned_bam"
    )
    trim_dir <- list("fastp_log", "unpaired")
    diff_dir <- list("DESeq2", "edgeR")
    pathway_types <- list("KEGG", "GO")
    kegg_types <- list("signalling", "metabolism", "disease", "sig_n_met")
    go_types <- list(
        "biological_process", "molecular_function", "cellular_component" )
    lapply(folder_to_create, checkcretdir, parentname = result.dir)
    lapply(trim_dir, checkcretdir, parentname = file.path(
        result.dir, "fastp_results"))
    lapply(diff_dir, checkcretdir, parentname = file.path(
        result.dir, "differential_analysis"))
    lapply(pathway_types, checkcretdir, parentname = file.path(
        result.dir, "gage_results"))
    lapply(kegg_types, checkcretdir, parentname = file.path(
        result.dir, "gage_results", "KEGG"))
    lapply(go_types, checkcretdir, parentname = file.path(
        result.dir, "gage_results", "GO"))
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
    return(c(
        qc.dir, trim.dir, deseq2.dir,
        edger.dir, gage.dir
    ))
}
    
