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
#' @param gsets : gene sets to analyse
#' @param work.dir : directory where results will be stored
#' @param same.dir : if the direction is same for GAGE analysis, GAGE parameter
#' @param compare : GAGE parameter
#' @param entity scientific name of organism
#' @param cpd_data compound data to run compound set analysis when available
#' @param gene_data gene data to run gene set analysis when available
#' @param cpd.idtype what type of id is cpd_data in 
#' @param gene_id_type what type of id is gene_id_type in
#' @import gage
#' @import utils
#' @import gage
#' @import pathview
#' @return nothing returned
#'


run_gage <- function(gsets, work.dir, same.dir,compare ,cpd_data, entity,
                    gene_data ,cpd.idtype = "kegg", gene_id_type ="entrez") {
#referenced from codes for pathview web available at https://pathview.uncc.edu/
    qcut<-0.01
    gage.dir <- dirname(dirname(work.dir))
    if (!is.null(cpd_data )){ fc_matrix <- cpd_data
        datatype <- "cpd"
    } else if (!is.null(gene_data)){ fc_matrix <- gene_data
        datatype <- "gene"
    } else if(!is.null(cpd_data) & !is.null(gene_data )){ datatype <- "combined"
    } else{ message("Please make sure one of the data type is present")}
    fc.kegg.p <- gage( fc_matrix, gsets = gsets, ref = NULL,samp = NULL,
                    same.dir = same.dir, compare = compare)
    pgs<-fc.kegg.p$greater[,"p.val"]
    pms<-cbind(pgs*2, (1-pgs)*2)
    pgs.gene<-apply(pms, 1, function(x) min(x))
    qgs.gene<-p.adjust(pgs.gene, method = "BH")
    colnames(pms)<-c("p.up", "p.dn")
    gage.out<-cbind(fc.kegg.p$greater[, c(2,5)], pms/2, p.val=pgs.gene, 
                    q.val=qgs.gene)
    gage.out<- gage.out[order(pgs.gene, -gage.out[,"set.size"]),]
    write.table(gage.out, file = file.path(work.dir , paste0("gage.out.", 
    datatype , ".tsv")), sep="\t", col.names=NA, quote = FALSE)
    sig.i<- gage.out[,"q.val"]<qcut & !is.na(gage.out[,"q.val"])
    nsig<-sum(sig.i, na.rm=TRUE)
    if(nsig>0) {
        gage.out.sig<-data.frame(gage.out)[sig.i,]
        ord1<-order(gage.out.sig[,"stat.mean"], decreasing=TRUE)
        gage.out.sig<-gage.out.sig[ord1,]
        write.table(gage.out.sig, file =file.path(work.dir, 
        paste0("gage.out.sig.", datatype , ".tsv")), sep="\t", 
                    col.names=NA, quote = FALSE)
        gpath_ids<-rownames(gage.out.sig)
    } else {
        message(paste0("No ", datatype,  " set selected in GAGE test, top 3 ", 
                    datatype,  "set plotted instead!", collapse = ""))
    gsel<-gage.out[,"set.size"]>0
    if(sum(gsel)>0) gpath_ids<-rownames(gage.out)[seq_len(min(sum(gsel),3))] }
    gpath_ids <- substr(gpath_ids, start = 4, stop = 8)
    if (same.dir == FALSE){
        plotpathways(gage.dir = work.dir, entity, gpath_ids[seq_len(6)] , 
                    cpd_data = cpd_data, gene_data =gene_data ,
                    cpd.idtype = "kegg", gene_id_type ="entrez")}
        return(list("pathways_selected"= gpath_ids,"pgs.gene" = pgs.gene,
                "gage_result"= gage.out,"data_used" =fc_matrix , 
                "gene_sets"=gsets ))
}
