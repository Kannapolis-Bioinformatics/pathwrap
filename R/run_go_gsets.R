#' run go enrichmnet analysis
#' 
#' This function runs go enrichment analysis for biological process, cellular component and metabolic function separately
#' 
#' @param gene_data logfold change, names should be genes
#' @param cnts gene counts 
#' @param outdir main directory for storing output of the process
#' @param entity Scientific name of species whose RNA is being analyzed
#' @param compare how the comparision is made in gage
#' @param use.fold if fold change should be used instead
#' @import pathview
#' @import gage
#' @param phenofile  file where the path of raw data is stored, used to generate ref and samp index



run_go_set_analysis <- function(gene_data,cnts, outdir,entity,compare,use.fold, phenofile ){
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
      gage.out  <- run_gage2(gene_data = gene_data, go_gene_sets[[i]], work.dir = outdir_list[i], same.dir = TRUE,
                                   compare = compare,gene_id_type = "ENTREZ", ref=NULL, samp=NULL)
    
    
    #################
    #TO DO
    #March 5
    #1) make the plot, find ref and samp index from phenofile
    #2) run compound
    #3) run combined compound and gene (code + thing)
    #run this later
      if (FALSE){ 
        phenofile_res <- process_phenofile(phenofile)
        coldata <- phenofile_res$coldata
        SampleName <- phenofile_res$SampleName
        filenames <- phenofile_res$FileName
        paired_info <- as.factor(as.character(phenofile_res$paired_info))
        
        cnts <- cnts[, coldata$SampleName]
        if (all(coldata$SampleName == colnames(cnts))) { # if this then proceed
          ref <- which(coldata$Class == levels(as.factor(coldata$Class))[1])
          samp <- which(coldata$Class == levels(as.factor(coldata$Class))[2])
        } else {
          message("make sure pheno file have only samples analysed")
        }
        
        gs=unique(unlist(kegg.gs[rownames(gse16873.kegg.p$greater)[1:3]]))
        
      for (gs in rownames(gage.out$less)[1:3]) {
           outname = gsub(" |:|/", "_", substr(gs, 12, 100))   #outdir_list[i]
      geneData(genes = go.sets.hs[[gs]], exprs = cnts, ref = ref, 
               samp = samp, outname = outname, txt = T, heatmap = T, limit = 3, 
               scatterplot = T)
      }
    }
}

    #fc_matrix, gsets, ref, samp, same.dir, compare, use.fold , work.dir, gcompare
}