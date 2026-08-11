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
    
    go_gene_sets <- list("biological_process" = go.bp, "cellular_component"=go.cc, "molecular_function" = go.mf)
    outdir_list <- list.files(file.path(outdir, "gage_results", "GO"), full.names = TRUE)
    
    #run gage for different go ontologies
    for (i in 1:length(go_gene_sets )){
      message("running for ", names(go_gene_sets)[i] )
      gage.out  <- run_gage2(gene_data = gene_data$logfoldchange, go_gene_sets[[i]], 
                              work.dir = outdir_list[i], same.dir = TRUE,
                                   compare = compare,gene_id_type = "ENTREZ", 
                              ref=NULL, samp=NULL,q_cutoff)
      
    
      plot_genedata( gene_data = gene_data$normalized_count, gage.out, 
                     gset =go_gene_sets[[i]] , outdir = outdir_list[i], 
                     mref=mref, msamp = msamp, compare= compare,q_cutoff )
      }
}

    #fc_matrix, gsets, ref, samp, same.dir, compare, use.fold , work.dir, gcompare
