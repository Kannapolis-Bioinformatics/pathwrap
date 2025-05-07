#' run KEGG enrichmnet analysis
#' 
#' This function runs go enrichment analysis for KEGG disease, metabolism and signalling. 
#' An additional analysis is also done by combining kegg metabolism and signalling pathways
#' 
#' @param gene_data dataframe or matrix or list of genes
#' @param outdir directory to store output
#' @param entity organism of interest
#' @param compare how the comparision is made
#' @import gage

run_kegg_analysis <- function(gene_data,outdir, entity, compare){
    ##download gene set 
        kegg.gs <- kegg.gsets(entity, check.new = TRUE)
        signmetinkegg <- kegg.gs$kg.sets[kegg.gs$sigmet.idx]
        diseaseinkegg <- kegg.gs$kg.sets[kegg.gs$dise.idx]
        siginkegg <- kegg.gs$kg.sets[kegg.gs$sig.idx]
        metainkegg <- kegg.gs$kg.sets[kegg.gs$met.idx]
        gene_sets = list( diseaseinkegg, metainkegg,signmetinkegg,siginkegg)
    
    ###run gage for different genesets
    outdir_list <- list.files(file.path(outdir, "gage_results", "KEGG"), full.names = T)
    for (i in 1:length(gene_sets)) {
        #run gage
        gpath_ids <- run_gage2( gene_data = gene_data, gene_sets[[i]], work.dir = outdir_list[i], same.dir = FALSE,
                compare = compare,gene_id_type = "ENTREZ",  ref=NULL, samp=NULL)
        print(gpath_ids["pids"])
        #plot sig pathways
        plotpathways(kegg.dir = outdir_list[i], entity, gpath_ids$pids[1:6] , 
                 cpd_data = NULL, gene_data =gene_data ,
                  gene_id_type ="entrez")
        
    }
    message("3c Mapping and plotting gene data to KEGG pathways complete." )
}
