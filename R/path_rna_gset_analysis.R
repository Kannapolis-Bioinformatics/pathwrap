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
#'
#' @param entity : scientific name of the species
#' @param fc_matrix : log fold expression values
#' @param compare : how the comparison is done for GAGE, see gage for details
#' @param outdir : directory in which GAGE will be run
#' @param cnts : counts of genes to use in pathview
#' @param entity organism of interest, scientific name
#' @param phenofile : path to the phenofile where raw data files path is stored 
#'
#' @import stats
#' @import utils
#' @import gage
#' @import pathview
#' @return gage_return_obj gage_return_obj which is list is returned
#'
#' @export

run_gene_gsets_analysis <- function(fc_matrix,cnts,outdir, entity, compare, 
                                         phenofile  ) {   
    

    #---------------GO ANALYSIS-----------------------------------------------#
    ####
    gage_go_res_file <- file.path(outdir, "gage_results", "GO", 
                                  "biological_process", "gage.out_greater.sig.gene.tsv")
    if (!file.exists(gage_go_res_file)){ #can be better
        run_go_set_analysis(fc_matrix,cnts, outdir, entity,  compare, phenofile) 
    }
    
    message("GO results are in ", file.path(outdir, "gage_results", "GO"))
    message("STEP 3a. GO ontology based GAGE analysis complete.")

    #----------------KEGG analysis---------------------------------------------#
    
    run_kegg_analysis (fc_matrix,outdir, entity,  compare )
    message("KEGG results are in ", file.path(outdir, "gage_results", "KEGG"))
    message("STEP 3b. KEGG pathway based GAGE analysis complete.")
    
  
}
