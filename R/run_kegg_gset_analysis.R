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
    ##download gene set 
        kegg.gs <- kegg.gsets(entity, check.new = TRUE)
        signmetinkegg <- kegg.gs$kg.sets[kegg.gs$sigmet.idx]
        diseaseinkegg <- kegg.gs$kg.sets[kegg.gs$dise.idx]
        siginkegg <- kegg.gs$kg.sets[kegg.gs$sig.idx]
        metainkegg <- kegg.gs$kg.sets[kegg.gs$met.idx]
        gene_sets = list( diseaseinkegg, metainkegg,signmetinkegg,siginkegg) #this is arranged/ordered
    
    ###run gage for different genesets
    outdir_list <- list.files(file.path(outdir, "gage_results", "KEGG"), full.names = T)
    for (i in 1:length(gene_sets)) {
        #run gage
        #names(gene_data$logfoldchange)<- rownames(gene_data$logfoldchange)
        #gpath_ids <- run_gage2( gene_data = gene_data$logfoldchange, gene_sets[[i]], work.dir = outdir_list[i], same.dir = FALSE,
        #        compare = compare,gene_id_type = "ENTREZ",  ref=NULL, samp=NULL)
        #print(gpath_ids["pids"])
        #plot sig pathways
        
        
        
        gpath_ids <- run_gage2( gene_data = logfoldchange, gene_sets[[i]], work.dir = outdir_list[i], same.dir = FALSE,
                                compare = compare,gene_id_type = "ENTREZ",  ref=NULL, samp=NULL, q_cutoff=q_cutoff)
        
        
        plotpathways(kegg.dir = outdir_list[i], entity, gpath_ids$pids[1:6] , 
                     cpd_data = NULL, gene_data =logfoldchange ,
                     gene_id_type ="entrez")
     
        
        print("this is ref col")
        print(mref)
        print("this is sam col")
        print(msamp)
        
        
       plot_genedata( gene_data = normalized_cnts, gpath_ids, 
                                     gset =gene_sets[[i]], outdir = outdir_list[i], 
                                     mref=mref, msamp = msamp, compare= compare,q_cutoff )
        
    }
    message("3c Mapping and plotting gene data to KEGG pathways complete." )
}
