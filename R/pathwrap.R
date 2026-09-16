#' Wrapper for RNASeq data analysis from raw reads to pathway visualization
#'
#' Wrapper that does quality control analysis of raw files, performs adapter and
#' quality trimming, builds genome index and does alignment, counts genes,
#' performs differential gene expression analysis, does gene enrichment testing using
#' GAGE, and visualizes enriched pathways using pathview, all in one
#' wrapper function. It can resume analysis if halted at any stage, generate
#' quality pictures, and produce comprehensive analysis of the data.
#'
#' @param phenofile file where the path of raw data is stored, its class and pair labels
#' @param entity Scientific name of species ("Homo sapiens" or "Mus musculus")
#' @param corenum number of cores available for run
#' @param ref.dir path to reference directory (*.fa and *.gtf); can be NULL or NA
#' @param cacheDir directory where temporary files created during alignment are stored
#' @param outdir main directory for storing output of the process
#' @param startover force complete rerun from beginning, default FALSE
#' @param aligner One of "Rhisat2" or "Rbowtie2"; Rbowtie2 can be very slow for
#'   eukaryotic species
#' @param diff.tool Differential expression tool to use: "DESeq2" or "edgeR"
#' @param cnts optional file path or pre-computed count matrix to skip alignment
#' @param gcompare comparison design for genes ("unpaired" or "paired")
#' @param npca number of top significant genes (by adjusted p-value) used for PCA
#' @param nheatmap number of genes for heatmap
#' @param fc_matrix logfold change matrix, used in place of running differential
#'   expression analysis when supplied (not currently implemented)
#' @param cdatapath data path for compound/metabolite data
#' @param csamp index/row numbers for samples, e.g. c(4, 5, 6)
#' @param cref index/row numbers for references, e.g. c(1, 2, 3)
#' @param ccompare compound comparison design ("unpaired" or "paired")
#' @param pathids pathway(s) of interest besides those enriched, e.g. c("04110"),
#'   c("04110", "04120"); set to NA to skip
#' @param nchunks number of chunks to run alignment, default 1
#' @param keep_tmp whether to keep trimmed fastq and intermediate BAM files, default TRUE
#' @param cpd.idtype compound id type (e.g., "KEGG COMPOUND accession", "KEGG")
#' @param q_cutoff cutoff for q-value used in pathway selection
#'
#' @return (invisibly) a list with elements:
#'   \code{status}, \code{outdir}, \code{gene_data}, \code{gene_path_results},
#'   \code{compound_data}, \code{manifest}
#' @export
pathwrap <- function(phenofile, entity, corenum = parallel::detectCores(),
                     ref.dir = NULL, cacheDir = tempdir(), outdir = NULL, startover = FALSE,
                     aligner = "Rhisat2", diff.tool = "DESeq2", cnts = NA, gcompare = "unpaired",
                     npca = 19, nheatmap = 10, fc_matrix = NA,
                     cdatapath = NA, csamp = NULL,
                     cref = NULL, ccompare = "paired", pathids = NA,
                     nchunks = 1, keep_tmp = TRUE, cpd.idtype = "KEGG", q_cutoff = 0.1) {
    
    on.exit(closeAllConnections(), add = TRUE)
    
    # -------------------------------------------------------------------------
    # 1. System & Input Prerequisites
    # -------------------------------------------------------------------------
    if (getRversion() < "4.4.0") {
        stop(
            "pathwrap requires R >= 4.4.0. You are running R ", as.character(getRversion()), ".\n",
            "Key dependencies (Rhisat2, txdbmaker) are unavailable in older R versions."
        )
    }
    
    if (!entity %in% c("Homo sapiens", "Mus musculus")) {
        stop("Parameter 'entity' must be either 'Homo sapiens' or 'Mus musculus'.")
    }
    
    if (!file.exists(phenofile)) {
        stop("Phenotype file not found: ", phenofile)
    }
    
    # Will alignment actually run, or will it be bypassed by a supplied count matrix?
    cnts_supplied <- (is.data.frame(cnts) || is.matrix(cnts)) ||
        (is.character(cnts) && length(cnts) == 1 && !is.na(cnts) && file.exists(cnts))
    alignment_will_run <- !cnts_supplied
    
    # Will the compound/pathway integration step run?
    compound_step_will_run <- is.character(cdatapath) && length(cdatapath) == 1 &&
        !is.na(cdatapath) && file.exists(cdatapath)
    
    # Is a user-specified pathway plot request active?
    pathids_active <- !(length(pathids) == 1 && is.na(pathids))
    
    # Deterministic, filesystem-safe checkpoint key for a given pathids value.
    compute_pathids_key <- function(ids) { 
        key <- paste(sort(ids), collapse = "_")
        key <- gsub("[^A-Za-z0-9_]", "-", key)
        if (nchar(key) > 60) {
            key <- paste0("h", as.character(sum(utf8ToInt(key)) %% 1e8))
        }
        key
    }
    
    # -------------------------------------------------------------------------
    # 2. Stage Directories & Persistent Logger
    # -------------------------------------------------------------------------
    if (is.null(outdir)) outdir <- file.path(getwd(), "results")
    outdir <- createdir(outdir, startover)
    if (is.null(outdir)) return(invisible(NULL))
    
    ckpt_dir <- file.path(outdir, ".checkpoints")
    if (!dir.exists(ckpt_dir)) dir.create(ckpt_dir, recursive = TRUE)
    if (!dir.exists(cacheDir)) dir.create(cacheDir, recursive = TRUE)
    
    log_file <- file.path(outdir, "pathwrap.log")
    write_log <- function(msg, console = TRUE) {
        timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
        entry <- paste(timestamp, msg)
        if (console) message(msg)
        cat(entry, file = log_file, sep = "\n", append = TRUE)
    }
    
    # -------------------------------------------------------------------------
    # 3. Stage Parameters & Detect Run Parity
    # -------------------------------------------------------------------------
    current_call <- as.list(environment())
    transient_vars <- c("startover", "write_log", "ckpt_dir", "log_file", "outdir",
                        "corenum", "pathids", "compute_pathids_key")
    
    # ref.dir/aligner/nchunks only matter if alignment actually runs
    if (!alignment_will_run) {
        transient_vars <- c(transient_vars, "ref.dir", "aligner", "nchunks")
    }
    # compound-related params only matter if the compound step actually runs
    if (!compound_step_will_run) {
        transient_vars <- c(transient_vars, "cpd.idtype", "csamp", "cref", "ccompare")
    }
    
    staged_params <- current_call[!names(current_call) %in% transient_vars]
    
    # Hash input files to detect data modifications
    files_to_hash <- c(phenofile = phenofile)
    if (is.character(cnts) && length(cnts) == 1 && !is.na(cnts) && file.exists(cnts)) {
        files_to_hash["cnts"] <- cnts
    }
    if (compound_step_will_run) {
        files_to_hash["cdatapath"] <- cdatapath
    }
    input_hashes <- vapply(files_to_hash, tools::md5sum, FUN.VALUE = character(1), USE.NAMES = TRUE)
    
    manifest_path <- file.path(outdir, ".run_manifest.rds")
    
    if (file.exists(manifest_path) && !startover) {
        prev_manifest <- readRDS(manifest_path)
        params_match <- isTRUE(all.equal(staged_params, prev_manifest$params))
        files_match  <- identical(input_hashes, prev_manifest$input_hashes)
        
        if (params_match && files_match) {
            write_log("Notice: Parameters and raw inputs match previous run.")
            
            final_cache <- file.path(ckpt_dir, "final_results.rds")
            
            # The full-cache shortcut is only valid if there isn't an active,
            # not-yet-rendered pathids request. pathids is intentionally excluded
            # from staged_params, so it is checked separately here.
            step4_pending <- FALSE
            if (pathids_active) {
                step4_file <- file.path(ckpt_dir, paste0("4_user_pathway_plot_", compute_pathids_key(pathids), ".rds"))
                step4_pending <- !file.exists(step4_file)
            }
            
            if (file.exists(final_cache) && !step4_pending) {
                write_log("Pipeline previously finished in full. Returning cached results.")
                return(invisible(readRDS(final_cache)))
            }
            write_log("Resuming partially completed pipeline from checkpoints...")
        } else {
            mismatches <- c()
            
            if (!params_match) {
                old_p <- prev_manifest$params
                new_p <- staged_params
                
                changed_names <- intersect(names(old_p), names(new_p))
                changed <- Filter(function(nm) !isTRUE(all.equal(old_p[[nm]], new_p[[nm]])), changed_names)
                
                added   <- setdiff(names(new_p), names(old_p))
                removed <- setdiff(names(old_p), names(new_p))
                
                detail <- c()
                if (length(changed) > 0) {
                    for (nm in changed) {
                        old_val <- paste(format(old_p[[nm]]), collapse = ", ")
                        new_val <- paste(format(new_p[[nm]]), collapse = ", ")
                        detail <- c(detail, paste0(nm, ": '", old_val, "' -> '", new_val, "'"))
                    }
                }
                if (length(added) > 0) {
                    detail <- c(detail, paste0("newly staged (branch now active): ", paste(added, collapse = ", ")))
                }
                if (length(removed) > 0) {
                    detail <- c(detail, paste0("no longer staged (branch now inactive): ", paste(removed, collapse = ", ")))
                }
                
                mismatches <- c(mismatches, paste0("parameter settings changed [", paste(detail, collapse = "; "), "]"))
            }
            
            if (!files_match) {
                old_h <- prev_manifest$input_hashes
                new_h <- input_hashes
                
                changed_files <- intersect(names(old_h), names(new_h))
                changed_files <- Filter(function(nm) !identical(old_h[[nm]], new_h[[nm]]), changed_files)
                
                added_files   <- setdiff(names(new_h), names(old_h))
                removed_files <- setdiff(names(old_h), names(new_h))
                
                file_detail <- c()
                if (length(changed_files) > 0) {
                    file_detail <- c(file_detail, paste0("content changed: ", paste(changed_files, collapse = ", ")))
                }
                if (length(added_files) > 0) {
                    file_detail <- c(file_detail, paste0("newly supplied: ", paste(added_files, collapse = ", ")))
                }
                if (length(removed_files) > 0) {
                    file_detail <- c(file_detail, paste0("no longer supplied: ", paste(removed_files, collapse = ", ")))
                }
                
                mismatches <- c(mismatches, paste0("input file contents [", paste(file_detail, collapse = "; "), "]"))
            }
            
            write_log(paste("Run inputs modified:", paste(mismatches, collapse = " | "), "- invalidating checkpoints."))
            unlink(list.files(ckpt_dir, full.names = TRUE))
        }
    } else if (startover) {
        write_log("Parameter 'startover=TRUE' passed. Clearing checkpoints.")
        unlink(list.files(ckpt_dir, full.names = TRUE))
    }
    
    # Save/update active run manifest
    saveRDS(list(params = staged_params, input_hashes = input_hashes, start_time = Sys.time()), manifest_path)
    
    # Checkpoint helper closure
    get_or_run <- function(stage_name, runner_fn) {
        stage_file <- file.path(ckpt_dir, paste0(stage_name, ".rds"))
        if (!startover && file.exists(stage_file)) {
            write_log(paste("Skipping completed step (cached):", stage_name))
            return(readRDS(stage_file))
        }
        result <- runner_fn()
        saveRDS(result, stage_file)
        return(result)
    }
    
    # -------------------------------------------------------------------------
    # 4. Processing Steps
    # -------------------------------------------------------------------------
    phenofile_res <- process_phenofile(phenofile)
    
    # Step 1a: FASTQC
    res.fastqc <- get_or_run("1a_fastqc", function() {
        res <- run_qc(phenofile_res$fq.dir, outdir, corenum)
        if (is.null(res)) {
            stop("One or more samples failed all QC metrics. Review raw fastq files.")
        }
        write_log("STEP 1a: FASTQC complete")
        return(res)
    })
    
    # Step 1b: FASTP
    get_or_run("1b_fastp", function() {
        run_fastp(phenofile, outdir, corenum)
        write_log("STEP 1b: FASTP complete")
        return(TRUE)
    })
    
    # Step 1c: Alignment & Counting (Skipped if cnts provided)
    aligned_obj_n_cnts <- get_or_run("1c_alignment_counting", function() {
        if (cnts_supplied) {
            write_log("Bypassing alignment: Using user-supplied count matrix.")
            counts_data <- if (is.character(cnts)) read.table(cnts, header = TRUE, row.names = 1) else cnts
            return(list(cnts = counts_data))
        }
        
        aln_res <- do_alignment_ncounting(phenofile_res, ref.dir, outdir,
                                          entity, cacheDir, aligner, corenum, nchunks)
        if (is.null(aln_res)) {
            stop("Alignment failed. Please check your reference files and indices.")
        }
        write_log("STEP 1c: ALIGNMENT AND COUNTING complete")
        return(aln_res)
    })
    
    # Step 2a: Differential Gene Expression
    gene_data <- get_or_run("2a_differential_expression", function() {
        res <- run_differerntial_gene_analysis(
            cnts = as.data.frame(aligned_obj_n_cnts$cnts),
            phenofile_res = phenofile_res,
            outdir = outdir,
            entity = entity,
            gcompare = gcompare,
            npca = npca,
            nheatmap = nheatmap,
            diff.tool = diff.tool
        )
        if (is.null(res)) stop("Differential gene expression analysis failed.")
        write_log(paste("STEP 2a: DIFFERENTIAL ANALYSIS complete using", diff.tool))
        return(res)
    })
    
    if (!keep_tmp) {
        delete_tmp_files(file.path(outdir, "fastp_results"))
        delete_tmp_files(file.path(outdir, "aligned_bam"))
    }
    
    # Step 2b: Gene Set Enrichment (GAGE)
    allpathresults <- get_or_run("2b_gene_gsets", function() {
        res <- run_gene_gsets_analysis(
            gene_data = gene_data,
            outdir = outdir,
            entity = entity,
            compare = gcompare,
            phenofile_res = phenofile_res,
            q_cutoff = q_cutoff
        )
        write_log("STEP 2b: GAGE enrichment complete")
        return(res)
    })
    
    # Step 3: Compound Integration & Pathview Visualization (Optional)
    if (compound_step_will_run) {
        # Step 3a: Compound KEGG gene set analysis
        cpd_res <- get_or_run("3a_compound_kegg_gsets", function() {
            res <- run_compound_kegg_gsets(
                cdatapath = cdatapath,
                cpd_id_type = cpd.idtype,
                csamp = csamp,
                keggorgcode = unname(gene_data$keggorgcode),
                cref = cref,
                ccompare = ccompare,
                outdir = outdir,
                q_cutoff = q_cutoff
            )
            write_log("STEP 3a: Compound KEGG gene sets analysis complete.")
            return(res)
        })
        
        # Step 3b: Combined pathway analysis (genes + compounds)
        path_hits <- get_or_run("3b_combined_pathway_analysis", function() {
            res <- run_combinedpath_analysis(
                outdir = outdir,
                gene_data = gene_data$logfoldchange,
                entity = entity,
                gcompare = gcompare,
                gage.out.cpd_res = cpd_res$gage.out.cpd,
                q_cutoff = q_cutoff
            )
            write_log("STEP 3b: Combined pathway analysis complete.")
            return(res)
        })
        
        # Step 3c: Pathview rendering
        # cpd_id_type is hardcoded to "kegg" here because cpd_res$cpd_data was
        # already converted to KEGG accessions inside run_compound_kegg_gsets().
        get_or_run("3c_pathview_rendering", function() {
            cpd_toplot <- as.matrix(cpd_res$cpd_data, rownames.force = TRUE)
            
            plotpathways(
                kegg.dir = file.path(outdir, "gage_results", "combined_analysis_kegg"),
                entity = entity,
                path.ids = path_hits,
                gene_data = gene_data$logfoldchange, cpd_data = cpd_toplot,
                gene_id_type = "entrez", cpd_id_type = "kegg")
            
            write_log("STEP 3c: Pathview rendering complete.")
            return(TRUE)
        })
        
        gage.out.cpd_res <- list(cpd_res = cpd_res, pathids = path_hits)
    } else {
        gage.out.cpd_res <- NULL
    }
    
    # -------------------------------------------------------------------------
    # Step 4: User-specified pathway plot(s) (optional)
    # -------------------------------------------------------------------------
    if (pathids_active) {
        pathids_key <- compute_pathids_key(pathids)
        
        get_or_run(paste0("4_user_pathway_plot_", pathids_key), function() {
            write_log(paste("Plotting user-specified pathway(s):", paste(pathids, collapse = ", ")))
            
            cpd_toplot   <- NULL
            cpd_id_type_ <- NULL
            if (compound_step_will_run) {
                cpd_toplot   <- as.matrix(cpd_res$cpd_data, rownames.force = TRUE)
                cpd_id_type_ <- "kegg"
            }
            
            plotpathways(
                kegg.dir     = file.path(outdir, "gage_results", "user_pathways", pathids_key),
                entity       = entity,
                path.ids     = pathids,
                gene_data    = gene_data$logfoldchange,
                cpd_data     = cpd_toplot,
                gene_id_type = "entrez",
                cpd_id_type  = cpd_id_type_
            )
            
            write_log("STEP 4: User-specified pathway plot(s) complete.")
            return(TRUE)
        })
    }
    
    final_output <- list(
        status = "completed",
        outdir = outdir,
        gene_data = gene_data,
        gene_path_results = allpathresults,
        compound_data = gage.out.cpd_res,
        manifest = readRDS(manifest_path)
    )
    
    saveRDS(final_output, file.path(ckpt_dir, "final_results.rds"))
    write_log("pathwrap pipeline execution completed successfully.")
    return(invisible(final_output))
}