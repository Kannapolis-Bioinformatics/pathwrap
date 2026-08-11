#' Run GAGE and Pathview
#'
#' `run_pathway` runs GAGE for pathway analysis
#' GAGE is based upon the user supplied comparison
#' method for the species specified.
#' The biological process, cellular component and
#' molecular function analysis for GO terms are done separately.
#'
#' KEGG disease, KEGG signalling and metabolism pathways are analysed separately
#' Top enriched pathways with  "q.val" < 0.1 are visualized using pathview.
#'
#' @param entity : scientific name of the species
#' @param gene_data : gene data with gene_data$logfoldchange #log fold expression values
#' @param compare : how the comparison is done for GAGE, see gage for details
#' @param outdir : directory in which GAGE will be run
#' @param entity organism of interest, scientific name
#' @param phenofile : path to the phenofile where raw data files path is stored 
#' @param q_cutoff q val for pathway selection
#' @import stats
#' @import utils
#' @import gage
#' @import pathview
#' @return gage_return_obj gage_return_obj which is list is returned
#'
#' @export

run_gene_gsets_analysis <- function(gene_data, outdir, entity, compare, 
                                         phenofile,q_cutoff  ) {   
    
  fc_matrix <- gene_data$logfoldchange
    ####
  plot.genedata = TRUE
  if (plot.genedata == TRUE) {
    phenofile_res <- process_phenofile(phenofile)
    coldata <- phenofile_res$coldata
    sampleName <- phenofile_res$sampleName
    filenames <- phenofile_res$FileName
    paired_info <- as.factor(as.character(phenofile_res$paired_info))
    if (all(coldata$sampleName == colnames(gene_data$normalized_count))) { # if this then proceed
      mref <- which(coldata$Class == levels(as.factor(coldata$Class))[1])
      msamp <- which(coldata$Class == levels(as.factor(coldata$Class))[2])
    } else {
      message("make sure pheno file have only msamples analysed")
    }
    
    
    gene_data$normalized_count <- gene_data$normalized_count[, coldata$SampleName]
    rownames(gene_data$normalized_count) <- str_remove(rownames(gene_data$normalized_count),pattern= "\\.\\d+")
    
    kegg.gs.species <- kegg.gsets(entity)
    if (sum(rownames(gene_data$normalized_count) %in% unlist(unname(kegg.gs.species$kg.sets))) < 10) {
      
      orgcode <- kegg.species.code(entity)
      org <- unname(korg[korg[,4]==entity, 3])
      twoletter = paste0(toupper(substr(org, 1, 1)), substr(org, 2, 2))
      entrezid <- id2eg( ids = rownames(gene_data$normalized_count),category="ENSEMBL",org=twoletter ) 
      gene_data$normalized_count <- mol.sum(gene_data$normalized_count, id.map = entrezid)
      #gene_data$normalized_count <- as.matrix(gene_data$normalized_count)
      names(gene_data$normalized_count) <- rownames(gene_data$normalized_count)
      
    }
  }
  #----------------------------------------------------------------------------#
  #---------------GO ANALYSIS-----------------------------------------------#
  #----------------------------------------------------------------------------#
    gage_go_res_file <- file.path(outdir, "gage_results", "GO", 
                                  "biological_process", "gage.out_greater.sig.gene.tsv")
    if (!file.exists(gage_go_res_file)){ #can be better
        run_go_set_analysis(gene_data , outdir, entity,  compare,use.fold=T, mref, msamp,q_cutoff) 
    }
    
    message("GO results are in ", file.path(outdir, "gage_results", "GO"))
    message("STEP 3a. GO ontology based GAGE analysis complete.")
  #----------------------------------------------------------------------------#
  #-----        KEGG analysis                         -----------------------#
  #----------------------------------------------------------------------------# 
  
  logfoldchange <- gene_data$logfoldchange
  normalized_cnts <-  gene_data$normalized_count
    run_kegg_analysis (logfoldchange, normalized_cnts ,  outdir, entity,  compare, mref, msamp,q_cutoff )

    message("KEGG results are in ", file.path(outdir, "gage_results", "KEGG"))
    message("STEP 3b. KEGG pathway based GAGE analysis complete.")
    
  
}
