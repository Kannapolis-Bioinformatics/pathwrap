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
#' @import gage
#' @import pathview
#' @return gage_return_obj_c return object from run_gage for compound


run_cpathway<- function( cdatapath,cpd_id_type= "KEGG COMPOUND accession",csamp,
                        cref, ccompare="paired" , cset_dir){
    if (grepl(".csv$", x = cdatapath)){
        cdata <- read.csv(cdatapath, header = T)
    } else{
        cdata<- read.table(file = cdatapath ,  sep="\t", header= TRUE)
    }    
    if (cpd_id_type!="KEGG COMPOUND accession"){
        cpd_idmap<-cpd2kegg(in.ids=rownames(cdata),in.type=toupper(cpd_id_type))
        didx<-duplicated(cpd_idmap[,1])
        cpd_idmap<-cpd_idmap[!didx,]

        cdata <-mol.sum(cdata, cpd_idmap)
    }
    
    # same.dir if to test for changes in a gene set toward a single direction
    if(!is.null(cref)){
        ncsamp<-length(csamp)
        ncref<-length(cref)
        if(ccompare=="paired" & ncsamp==ncref) 
            {cdata<-cdata[,csamp]- cdata[,cref]}
        else if (ncref==1) cdata<- cdata[,csamp]- cdata[,cref]
        else cdata <-cdata[,csamp]- rowMeans(cdata[,cref])
    }
    csets <- loadcsets()
    gage_return_obj_c <- run_gage(csets, cset_dir , same.dir = FALSE,
            compare = ccompare, fc_matrix = cdata,entity= NA)
    return(gage_return_obj_c)
    
    
}
