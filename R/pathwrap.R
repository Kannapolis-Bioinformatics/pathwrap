#' Wrapper for RNASeq data analysis from raw reads to pathway visualization
#'
#' wrapper that does quality control analysis of raw files, performs adapter and
#' quality trimming, builds genome index and does alignment , counts genes
#' perfoms differential gene expression analysis,does gene enrichment test using
#' GAGE and visualize the enriched pathways using pathview all using one
#' wrapper function. It has the ability to continue the analysis if it is
#' halted at any stage and generate quality pictures and generate comprehensive
#' analysis of the data.
#' 
#' @param outdir main directory for storing output of the process
#' @param startover do you want to start from beginning, default =False
#' @param phenofile file where the path of raw data is stored.
#' @param corenum number of cores avaialble for run
#' @param ref.dir path to reference directory which contain
#' reference file(*.fa) and annotation file(*.gtf), can be NA
#' @param entity Scientific name of species whose RNA is being analyzed
#' @param cacheDir directory where temporary files created during alignment
#' @param aligner One of "Rhisat2" or "Rbowtie2"; Rbowtie2 can be very slow
#'    for human and eukaryotic species
#' @param gcompare how the comparision is done for transcripts/genes
#' @param npca number of genes to use for pca
#' @param nheatmap number of genes for heatmap
#' @param fc_matrix logfold change for gage if deseq2 is not run
#' @param cdatapath data path for compound data
#' @param cpd_id_type "KEGG COMPOUND accession"
#' @param csamp index/row number where sample files are ex: c(5,6,4)
#' @param cref index/row number where references are, ex: c(1,2,3)
#' @param ccompare how the compound data is compared, default paired
#' @param qcut threshold for pathway selection, default 0.01
#' @param pathids pathway of interest only necessary if enrichment is not run
#' @param nchunks default 1, in how many chunks you want to run alignment
#' @param keep_tmp weather to store aligned bam files and trimmed fastq files default True
#' @param diff.tool weather to use "DESeq2" or edgeR for differential gene analysis
#' @param cnts file path to the cnts per gene from other analysis
#' @return returns message if analysis complete successfully
#' @export
pathwrap <- function(  phenofile,entity,corenum=detectCores(),ref.dir=NULL, cacheDir=tempdir(), 
                       outdir=NULL, startover=FALSE,
                     aligner="Rhisat2",diff.tool = "DESeq2", cnts=NA,   gcompare="unpaired", 
                     npca= 19, nheatmap=10,fc_matrix=NA, 
                     cdatapath=NA,cpd_id_type= "KEGG COMPOUND accession",csamp=NULL,
                     
                     cref=NULL, ccompare="unpaired" ,  qcut=0.01, pathids="04110",
                     nchunks=1, keep_tmp = TRUE,  cpd.idtype = "KEGG"){
    
    on.exit(closeAllConnections())
    #A. PREPARE DIRECTORIES
    
    
    outdir <- createdir(pos =1, outdir, startover)
   
    phenofile_res <- process_phenofile(phenofile)
    
    res.fastqc <- run_qc(phenofile_res$fq.dir, outdir, corenum)
    message("STEP 1a: FASTQC complete")
    if(is.null(res.fastqc)){
        message("One or more samples failed all the QC metrics.\n",
                "Please remove them from the raw directory and phenofile to continue.")
        return(invisible(NULL))
    }
    
    run_fastp(phenofile, outdir, corenum)
    message("STEP 1b: FASTP complete")
    aligned_obj_n_cnts <- do_alignment_ncounting(phenofile, ref.dir, outdir, entity, cacheDir, 
                                                 aligner, corenum,nchunks)
    message("STEP 1c: ALIGNMENT AND COUNTING complete" )
    print("this is cnts")
    print(dim(as.data.frame(aligned_obj_n_cnts$cnts)))
    #cnts can be filename of cnts
    #phenofile can contain only samplename and class
     gene_data <- run_differerntial_gene_analysis(cnts = as.data.frame(aligned_obj_n_cnts$cnts), phenofile, outdir, entity,  gcompare, npca, nheatmap, diff.tool = diff.tool)
     if(is.null(gene_data)){
         return(invisible(NULL))
     }
     message("STEP 2a: DIFFERENTIAL ANALYSIS complete using ", diff.tool)
     if (!keep_tmp){
         delete_tmp_files(file.path(outdir, "fastp_results"))
         delete_tmp_files(file.path(outdir, "aligned_bam"))
        
     }
    # #check from here
     run_gene_gsets_analysis(fc_matrix = gene_data$logfoldchange , as.data.frame(aligned_obj_n_cnts$cnts), outdir, entity,  gcompare,phenofile)#, use.fold=TRUE)
    
     if (!is.na(cdatapath)){
         gage.out.cpd_res <- run_compound_kegg_gsets(cdatapath,cpd_id_type= cpd.idtype,csamp,
                                                     cref, ccompare="paired" , outdir, entity)
        
         #determine how you want to combine 
         pathids<- run_combinedpath_analysis(outdir, gene_data$logfoldchange, entity, gcompare,gage.out.cpd_res$gage.out.cpd,qcut=qcut)
        
         plotpathways(kegg.dir = file.path(outdir, "gage_results", "combined_analysis_kegg"), entity, pathids , 
                      cpd_data = gage.out.cpd_res$cpd_data, gene_data =gene_data$logfoldchange ,
                      gene_id_type ="entrez", cpd_id_type=cpd.idtype)
     }
    return("the analysis complete successfully")
    
    
}
