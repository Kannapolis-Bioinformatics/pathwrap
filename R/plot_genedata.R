#' Plot the genes belonging to top significant pathways
#' whose gene expression is more than that of background genes
#' 
#' `plot_genedata` Plot the genes belonging to top significant pathways
#' whose expression is more than that of background genes for all the samples 
#' genes counts are normalized counts using deseq2
#' 
#' @param gene_data obj returned by differential analysis function
#' @param gage.out result of gage analysis for logFC data
#' @param gset gene set of interest
#' @param outdir directory to store the heatmaps
#' @param compare : GAGE parameter
#' @param mref reference index for gene/cpd data
#' @param msamp sample index for gene/cpd data
#' @param q_cutoff q val  to select pathways
#' @import gage
#' @import utils
#' @import gage
#' @import pathview
#' @return nothing returned
#'

plot_genedata <- function( gene_data, gage.out, gset, outdir, mref, msamp, compare,q_cutoff){
    fc.kegg.p <- gage.out$gage.out
    sel <- fc.kegg.p$greater[, "q.val"] < q_cutoff &  !is.na(fc.kegg.p$greater[, "q.val"])
    if(sum(sel)==0){
        sel <- c(1,2,3,4,5,6)
    }
    top_ids <- na.omit(rownames(fc.kegg.p$greater)[sel][1:6])
    gs      <- unique(unlist(gset[top_ids]))
    essData <- essGene(gs, gene_data, ref=mref, samp=msamp, compare=compare)
    #gs can be name of gene set , or gene id vector
    
    #essential member genes in a gene set.#change over noise
    #all genes in pathways  
    #for (id in na.omit(rownames(gage.out$gage.out$greater)[1:12])) {
    for (id in top_ids){
        genes_in_term=unique(unlist(gset[id]))
        outname <- file.path(outdir,  paste0(gsub(" |:|/", "_", id) , "_greater") )
        if (sum(genes_in_term %in% rownames(gene_data)) < 3) {
            message(paste("Skipping", id, "- not enough genes (less than 3) found."))
            next
        }
        geneData(genes = genes_in_term, exprs =essData , ref = mref, 
                 samp = msamp, outname = outname, txt = T, heatmap = T, limit = 3, 
                 scatterplot = T)
    }
    
    if (!is.null(gage.out$gage.out$less)){
        fc.kegg.p <- gage.out$gage.out
        sel <- fc.kegg.p$less[, "q.val"] < q_cutoff &  !is.na(fc.kegg.p$less[, "q.val"])
        if(sum(sel)==0){
            sel <- c(1,2,3,4,5,6)
        }
        top_ids <- na.omit(rownames(fc.kegg.p$less)[sel][1:6])
        gs      <- unique(unlist(gset[top_ids]))
        essData <- essGene(gs, gene_data, ref=mref, samp=msamp, compare=compare)
        #gs can be name of gene set , or gene id vector
        
        #essential member genes in a gene set.#change over noise
        #all genes in pathways  
        #for (id in na.omit(rownames(gage.out$gage.out$less)[1:12])) {
        for (id in top_ids){
            genes_in_term=unique(unlist(gset[id]))
            outname <- file.path(outdir, paste0( gsub(" |:|/", "_", id) , "_less") )
            if (sum(genes_in_term %in% rownames(gene_data)) < 3) {
                message(paste("Skipping", id, "- not enough genes (less than 3) found."))
                next
            }
            geneData(genes = genes_in_term, exprs =essData , ref = mref, 
                     samp = msamp, outname = outname, txt = T, heatmap = T, limit = 3, 
                     scatterplot = T)
        } 
    }
    
    
    
    ###
    non_redunant = FALSE
    if (non_redunant){
    
    check_and_warn <- function(condition, message) {
        if (condition) {
            warning(message)
        }
    }

            tryCatch({
                coreset <- esset.grp(gage.out$gage.out$greater, gene_data,gsetm ,  ref=mref, samp=msamp, compare=compare, use.q = T, cutoff = q_cutoff)
                for (id in na.omit(rownames(gage.out$gage.out$greater)[1:12])) {
                    genes_in_term=unique(unlist(gsetm[id]))
                    outname <- file.path(outdir,  paste0("core", gsub(" |:|/", "_", id) ) )
                    core_genes <- coreset$coreGeneSets[[id]]
                    
                    message("length of genes selected in pathway")
                    message(length(genes_in_term))
                    message(sum(genes_in_term %in% rownames(gene_data)))
                    message("above is common genes found")
                    if (sum(genes_in_term %in% rownames(gene_data)) < 3) {
                        message(paste("Skipping", id, "- not enough genes (less than 3) found."))
                        next
                    }
                    message("this is dim of essData")
                    message(dim(essData))
                    #exprs = gene_data [ names(gene_data)%in%coreset$coreGeneSets,]
                    geneData(genes = core_genes, exprs = gene_data, ref = mref, 
                             samp = msamp, outname = outname, txt = T, heatmap = T, limit = 3, 
                             scatterplot = T)
                    
                }
            }, error = function(w) {
                check_and_warn(TRUE,"there are no significant gene set, try lowering q val")
            })
    }
}
            
    
    
#         
#         #############################
#        coreset <- esset.grp(gage.out$gage.out$greater, gene_data,gsetm ,  ref=mref, samp=msamp, compare=compare, use.q = T, cutoff = q_cutoff)
#         for (id in na.omit(rownames(gage.out$gage.out$greater)[1:12])) {
#             genes_in_term=unique(unlist(gsetm[id]))
#             outname <- file.path(outdir,  paste0("core", gsub(" |:|/", "_", id) ) )
#             core_genes <- coreset$coreGeneSets[[id]]
#             
#             print("length of genes selected in pathway")
#             print(length(genes_in_term))
#             print(sum(genes_in_term %in% rownames(gene_data)))
#             print("above is common genes found")
#             if (sum(genes_in_term %in% rownames(gene_data)) < 3) {
#                 message(paste("Skipping", id, "- not enough genes (less than 3) found."))
#                 next
#             }
#             print("this is dim of essData")
#             print(dim(essData))
#             #exprs = gene_data [ names(gene_data)%in%coreset$coreGeneSets,]
#             geneData(genes = core_genes, exprs = gene_data, ref = mref, 
#                      samp = msamp, outname = outname, txt = T, heatmap = T, limit = 3, 
#                      scatterplot = T)
#             
#         }
#     
# }
    
