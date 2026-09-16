#' run go enrichmnet analysis
#' 
#' This function runs go enrichment analysis for biological process, cellular component and metabolic function separately
#' 
#' @param gene_data object returned from prev names should be genes
#' @param mref reference  colums for plotting gene data
#' @param msamp sample columns for plotting gene data
#' @param q_cutoff q value for pathway selection
#' @param outdir main directory for storing output of the process
#' @param entity Scientific name of species whose RNA is being analyzed
#' @param compare how the comparision is made in gage
#' @param use.fold if fold change should be used instead
#' @import pathview
#' @import gage




run_go_set_analysis <- function(gene_data, outdir,entity,compare,use.fold=T, mref, msamp,q_cutoff ){
    data(korg, package = "pathview")
    keggcode_sel <- unname(korg[which(korg[, 4] == entity), 3])
    data(bods, package = "gage", envir = environment())
    common_name_species <- bods[, 2][which(bods[, 3] == keggcode_sel)]
    
    #download go ontologies
    go.gs <- go.gsets(common_name_species)
    go.bp <- go.gs$go.sets[go.gs$go.subs$BP]
    go.mf <- go.gs$go.sets[go.gs$go.subs$MF]
    go.cc <- go.gs$go.sets[go.gs$go.subs$CC]
    
    go_gene_sets <- list("biological_process" = go.bp, "cellular_component"=go.cc, 
                         "molecular_function" = go.mf)
    go_folder_map <- c("biological_process" = "biological_process",
                       "cellular_component" = "cellular_component",
                       "molecular_function" = "molecular_function")
    go_base <- file.path(outdir, "gage_results", "GO")
    #go_all_results <- vector("list", length(go_gene_sets))
    #names(go_all_results) <- names(go_gene_sets)
    
    #run gage for different go ontologies
    for (category in names(go_gene_sets)) {
      message("running for ", category )
      work_dir <- file.path(go_base, go_folder_map[[category]])
      
      gage.out  <- run_gage2(gene_data = gene_data$logfoldchange, 
                             go_gene_sets[[category]], 
                              work.dir = work_dir ,same.dir = TRUE,
                                   compare = compare,gene_id_type = "ENTREZ", 
                              ref=NULL, samp=NULL,q_cutoff)
      
    
      plot_genedata( gene_data = gene_data$normalized_count, gage.out, 
                     gset =go_gene_sets[[category]] , outdir = work_dir ,
                     mref=mref, msamp = msamp, compare= compare,q_cutoff )
      #go_all_results[[category]] <- gage.out
      }
}

    #fc_matrix, gsets, ref, samp, same.dir, compare, use.fold , work.dir, gcompare
