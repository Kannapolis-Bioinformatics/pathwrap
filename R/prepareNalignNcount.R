#' Align and count
#' This function identifies the type of data(SE or PE), creates sampleFile neccessary for qAlign, 
#' checks if references are okay for alignement and counting, runs qAlign and qcount 
#' @param corenum number of cores to run alignment 
#' @param nchunks if the alignment should be done in chunks
#' @param phenofile : path to the phenofile where raw data files path is stored 
#' @param ref.dir directory where reference genome and annotation is stored
#' @param outdir directory for output files
#' @param entity organism of interest
#' @param cacheDir : directory where temporary files generated during alignment
#' @param aligner : weather Rhisat2 or Rbowtie should be used for alignment
#' @return lists with names cnts and aligned_proj 
#' @export
do_alignment_ncounting <- function(phenofile, ref.dir, outdir, entity, cacheDir, aligner, corenum,nchunks ){
    phenofile_res <- process_phenofile(phenofile)
    references <- check_references(ref.dir, outdir, entity)#, compare)
    if (is.null(references)){
        
        return(invisible(NULL))
    }
    message("References are Ok")
    aligned_proj_list <- run_qAlign(phenofile, cacheDir, aligner, references, outdir,corenum,nchunks)
    message("Alignment OK")
   qCountdf <- run_qCount(aligned_proj_list, corenum, outdir ,entity, references)
   if(! is.null(ref.dir)){
       onexistcleanup(ref.dir, entity)
   }
    return(list("cnts"= qCountdf))#, "aligned_proj"=aligned_proj ))
   #return(list("cnts"= qCountdf, "aligned_proj"=aligned_proj ))
}