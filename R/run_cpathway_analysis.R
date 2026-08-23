#' Function for running the compound pathway analysis
#'
#' This function prepares for running the gage analysis on compound sets
#' @param keggorgcode three char code for kegg organism
#' @param q_cutoff  q value cutoff for selecting pathway
#' @param cdatapath :data path for compound data
#' @param cpd_id_type type of compound id as compatible to ones in 
#' pathview::data(rn.list)
#' @param csamp column number of sample
#' @param cref column number of reference
#' @param ccompare comparison type for sample
#' @param outdir directory where results of compound set analysis is stored
#' @import gage
#' @import pathview
#' @return gage_return_obj_c return object from run_gage for compound
#' @export

run_compound_kegg_gsets<- function( cdatapath,cpd_id_type= "KEGG COMPOUND accession",csamp,keggorgcode,
                        cref, ccompare="paired" , outdir, q_cutoff){
#referenced from codes for pathview web  https://pathview.uncc.edu/
    cset_dir <- file.path(outdir, "gage_results", "kegg_csets")
    if (grepl(".csv$", x = cdatapath)){
        cpd_data <- read.csv(cdatapath, header = TRUE, row.names = 1)
    } else{
        cpd_data<- read.table(file = cdatapath,sep="\t",header= TRUE,
                            row.names = 1)
        }
    
    if (cpd_id_type =="KEGG COMPOUND accession"){cpd_id_type <- "KEGG" }
    if (cpd_id_type!="KEGG"){# COMPOUND accession" | cpd_id_type != "KEGG"){
        cpd_idmap<-cpd2kegg(in.ids=rownames(cpd_data),
                            in.type=toupper(cpd_id_type))
        didx<-duplicated(cpd_idmap[,1])
        cpd_idmap<-cpd_idmap[!didx,]
        cpd_data <-mol.sum(cpd_data, cpd_idmap)}
    
    # same.dir if to test for changes in a gene set toward a single direction
    if(!is.null(cref)){
        ncsamp<-length(csamp)
        ncref<-length(cref)
        if(ccompare=="paired" & ncsamp==ncref) 
            {cpd_data<-cpd_data[,csamp]- cpd_data[,cref]}
        else if (ncref==1) cpd_data<- cpd_data[,csamp]- cpd_data[,cref]
        else cpd_data <-cpd_data[,csamp]- rowMeans(cpd_data[,cref]) }
    csets <- loadcsets(keggorgcode)
    #names(csets)<-str_replace_all(names(csets),"^map",kegg.species.code(entity))
    message("gage running submitted for compound")
    gpath_ids <- run_gage2(gsets = csets,  same.dir = FALSE,
            compare = ccompare,gene_data= cpd_data, 
            , ref=NULL, samp=NULL, work.dir=cset_dir,q_cutoff=q_cutoff)
   
    
    
    
    plotpathways(kegg.dir =cset_dir, keggorgcode, gpath_ids$pids , 
                 cpd_data = cpd_data, gene_data = NULL ,
                 cpd_id_type= "kegg")
    
    return(list("pids"= gpath_ids$pids, "cpd_data"=cpd_data, "gage.out.cpd" = gpath_ids$gage.out))
}
