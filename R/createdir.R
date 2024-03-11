#' create the directories for saving the result
#'
#' @param entity : scientific name of the organism
#' @param startover : if the run is again?
#' @param keep_tmp : keep bam files resulting from alignment or not
#' @param pos enviromental 
#' @param outdir where the directory and its subfolder will be created
#' @return list of folder names
#' 
# qc.dir <- fastqc_results
# diff.dir <- differential_analysis
# trim.dir <- fastp_results
# gage.dir <- gage_results
# trim.log <- fastp_log
# edger.dir <- edgeR
# message(edger.dir)
# deseq2.dir <- DESeq2
# kegg.dir <- KEGG
# go.dir <- GO
#     

createdir <- function(pos =1, outdir, entity, startover, keep_tmp) {
    on.exit(closeAllConnections())
    if (file.exists(outdir)& length(list.files(outdir))>0 & startover == TRUE) {
        ans <- readline(paste0("Are you sure you want to delete everything in ",
                            outdir , "? "))
        if (substr(ans, 1, 1) == "n"){
            message("Make sure result directory is empty to start from beginning
                    or use startover = FALSE with same other parameters")
            return(invisible(x=NULL))}
        else {  unlink(outdir, recursive = TRUE)}}
    if (startover == FALSE){
        ans <- readline(paste0("Is this your first run and are parameters 
                            same as previous runs? "))
        if (substr(ans, 1, 1) == "n"){
            message("Make sure the parameters are same as previous run or answer
                    yes to above question ")
            return(invisible(x=NULL))}}
        if (!file.exists(outdir)) {
            dir.create(outdir)} # default output file
        message(paste0("The results will be organized in ", outdir, 
                    collapse = "|"))
        checkcretdir <- function(parentname, dirname) {
            if (!file.exists(file.path(parentname, dirname,
                                    fsep = .Platform$file.sep))) {
                dir.create(file.path(parentname, dirname,
                                    fsep = .Platform$file.sep))}
            assign(dirname,# .GlobalEnv)#environment())
                value = file.path(parentname, dirname,
                    fsep = .Platform$file.sep),envir = as.environment(pos))} 
        folder_to_create <- list("fastqc_results", "fastp_results",
                        "gage_results","differential_analysis", "aligned_bam")
        trim_dir <- list("fastp_log", "unpaired")
        diff_dir <- list("DESeq2", "edgeR")
        pathway_types <- list("KEGG", "GO", "KEGG_CSETS")
        kegg_types <- list("signalling", "metabolism", "disease", "sig_n_met")
        go_types <- list(
            "biological_process", "molecular_function", "cellular_component" )
        lapply(folder_to_create, checkcretdir, parentname = outdir)
        lapply(trim_dir, checkcretdir, parentname = file.path(
            outdir, "fastp_results",fsep = .Platform$file.sep))
        lapply(diff_dir, checkcretdir, parentname = file.path(
            outdir, "differential_analysis",fsep = .Platform$file.sep))
        lapply(pathway_types, checkcretdir, parentname = file.path(
            outdir, "gage_results",fsep = .Platform$file.sep))
        lapply(kegg_types, checkcretdir, parentname = file.path(
            outdir, "gage_results", "KEGG",fsep = .Platform$file.sep))
        lapply(go_types, checkcretdir, parentname = file.path(
            outdir, "gage_results", "GO",fsep = .Platform$file.sep))
        return(c(fastqc_results, fastp_results, DESeq2,edgeR, gage_results,
                KEGG_CSETS))
}
