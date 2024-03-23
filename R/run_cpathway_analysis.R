#' Function for running the compound pathway analysis
#'
#' This function prepares for running the gage analysis on compound sets
#'
#' @param cdatapath :pathway for compound data
#' @param cpd_id_type type of compound id as compatible to ones in 
#' pathview::data(rn.list)
#' @param csamp column number of sample
#' @param cref column number of reference
#' @param ccompare comparision type for sample
#' @param cset_dir directory where results of compound set analysis is stored
#' @param entity scientific name of the organism of study
#' @import gage
#' @import pathview
#' @return gage_return_obj_c return object from run_gage for compound


run_cpathway<- function( cdatapath,cpd_id_type= "KEGG COMPOUND accession",csamp,
                        cref, ccompare="paired" , cset_dir, entity){
#referenced from codes for pathview web available at https://pathview.uncc.edu/
    
    if (grepl(".csv$", x = cdatapath)){
        cpd_data <- read.csv(cdatapath, header = TRUE, row.names = 1)
    } else{
        cpd_data<- read.table(file = cdatapath,sep="\t",header= TRUE,
                            row.names = 1)}
    
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
    csets <- loadcsets()
    gage_return_obj_c <- run_gage(csets, cset_dir , same.dir = FALSE,
            compare = ccompare, entity= entity, cpd_data, 
            cpd.idtype = "kegg",gene_data=NULL, gene_id_type  =NA)
    return(gage_return_obj_c)
}
