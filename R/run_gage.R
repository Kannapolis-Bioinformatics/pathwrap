#' Run GAGE and Pathview
#'
#' `run_pathway` runs GAGE for pathway analysis
#' GAGE is based upon the user supplied comparision
#' method for the species specified.
#' The biological process, cellular component and
#' molecular function analysis for GO terms are done separately.
#'
#' KEGG disease, KEGG signalling and metabolism pathways are analysed separately
#' Top enriched pathways with  "q.val" < 0.1 are visualized using pathview.
#' @param q_cutoff cutoff of q value for pathway selection
#' @param gsets : gene sets to analyse
#' @param work.dir : directory where results will be stored
#' @param same.dir : if the direction is same for GAGE analysis, GAGE parameter
#' @param compare : GAGE parameter
#' @param gene_data gene data to run gene set analysis when available
#' @param gene_id_type what type of id is gene_id_type in
#' @param ref reference index for gene/cpd data
#' @param samp sample index for gene/cpd data
#' @import gage
#' @import utils
#' @import gage
#' @import pathview
#' @return nothing returned
#'

run_gage2 <- function(gene_data, gsets,  same.dir, compare, work.dir, gene_id_type = "ENTREZ", ref=NULL, samp=NULL, q_cutoff ){
    
    #message("referenced from codes for pathview web available at https://pathview.uncc.edu/
    #(same.dir = FALSE).  2d test (greater and stats)")
    
    #qcut <-0.01
    gage.dir <- dirname(dirname(work.dir))
    
    fc.kegg.p <- gage( gene_data, gsets = gsets, ref = ref,samp = samp,
                       same.dir = same.dir, compare = compare)
    message("gage ran successfully")
    write.table(fc.kegg.p$greater,file= file.path(work.dir, 
                                                  "gage.out_greater.sig.tsv"), sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
    
    #find top up regulated genes
    #matrix is sorted by global p value or q 
    sel <- fc.kegg.p$greater[, "q.val"] < q_cutoff &  !is.na(fc.kegg.p$greater[, "q.val"])
    message("the total number of "  , basename(work.dir), " pathways enriched is ", sum(sel))
    if(sum(sel)==0){
        sel <- c(1,2,3,4,5,6)
    }
    path.ids <-  as.character(na.omit(rownames(fc.kegg.p$greater)[sel][1:6]))
    
    if(same.dir == TRUE){
        write.table(fc.kegg.p$less,file= file.path(work.dir, 
                                                   "gage.out_less.sig.tsv"), 
                    sep = "\t", quote = F, col.names = NA, row.names = TRUE)
    
        sel <- fc.kegg.p$less[, "q.val"] < q_cutoff  &  !is.na(fc.kegg.p$less[, "q.val"])
        if(sum(sel)==0){
            sel <- c(1,2,3,4,5,6)
        }
        path.ids_less <-  as.character(na.omit(rownames(fc.kegg.p$less)[sel][1:6]))
        path.ids <- c(path.ids,path.ids_less )
        #return(path.ids) #GO:0140053 mitochondrial gene expression"                      
    }
    
    
    if (grepl(pattern = "GO",  x = work.dir)) { 
        path.ids2 <-  str_sub( path.ids , 4,10)
    }else {
        path.ids2 <-  str_sub( path.ids , 4,8)
        }
    
    message("this is path ids found by gage")
    message(path.ids2)
    return(list("pids"= path.ids2, "gage.out" = fc.kegg.p ))
    
    # return(list("pathways_selected"= gpath_ids,"pgs.gene" = pgs.gene,
    #         "gage_result"= gage.out,"data_used" =gene_data , 
    #         "gene_sets"=gsets ))
}
