#' run KEGG enrichmnet analysis
#' 
#' This function runs go enrichment analysis for KEGG disease, metabolism and signalling. 
#' An additional analysis is also done by combining kegg metabolism and signalling pathways
#' @param normalized_cnts normalized cnts to plot gene data 
#' @param mref reference  colums for plotting gene data
#' @param msamp sample columns for plotting gene data
#' @param q_cutoff q value for pathway selection
#' @param logfoldchange dataframe or matrix or list of genes
#' @param outdir directory to store output
#' @param entity organism of interest
#' @param compare how the comparision is made
#' @import gage

run_kegg_analysis <- function(logfoldchange, normalized_cnts , outdir, entity, compare, mref, msamp, q_cutoff){

        kegg.gs <- kegg.gsets(entity, check.new = TRUE)
        signmetinkegg <- kegg.gs$kg.sets[kegg.gs$sigmet.idx]
        diseaseinkegg <- kegg.gs$kg.sets[kegg.gs$dise.idx]
        siginkegg <- kegg.gs$kg.sets[kegg.gs$sig.idx]
        metainkegg <- kegg.gs$kg.sets[kegg.gs$met.idx]
        
        gene_sets = list( "disease" = diseaseinkegg, 
                          "met"=metainkegg,
                          "sigmet"=signmetinkegg,
                          "sig"=siginkegg) 
        

        folder_map <- c("disease" = "disease",
                        "met"     = "metabolism",
                        "sigmet"  = "signmet",
                        "sig"     = "signalling")
        kegg_base <- file.path(outdir, "gage_results", "KEGG")
        
        all_results <- vector("list", length(gene_sets))
        names(all_results) <- names(gene_sets)

    #outdir_list <- list.files(file.path(outdir, "gage_results", "KEGG"), full.names = T)
    for (category in names(gene_sets)) {
        work_dir <- file.path(kegg_base, folder_map[[category]])
        
        #gpath_ids <- run_gage2( gene_data = logfoldchange, gene_sets[[i]], work.dir = outdir_list[i], same.dir = FALSE,
        #                        compare = compare,gene_id_type = "ENTREZ",  ref=NULL, samp=NULL, q_cutoff=q_cutoff)
        gpath_ids <- run_gage2( gene_data = logfoldchange, gene_sets[[category]], work.dir = work_dir, same.dir = FALSE,
                                compare = compare,gene_id_type = "ENTREZ",  ref=NULL, samp=NULL, q_cutoff=q_cutoff)
        
        
        plotpathways(kegg.dir =work_dir, entity, gpath_ids$pids[1:6] , 
                     cpd_data = NULL, gene_data =logfoldchange ,
                     gene_id_type ="entrez")
        
       plot_genedata( gene_data = normalized_cnts, gpath_ids, 
                                     gset =gene_sets[[category]], outdir = work_dir, 
                                     mref=mref, msamp = msamp, compare= compare,q_cutoff )
       all_results[[category]] <- gpath_ids  
       
        
    }
    message("3c Mapping and plotting gene data to KEGG pathways complete." )
    return(all_results) #this is related to only one type of pathway
}
