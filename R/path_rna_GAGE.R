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
#' @param fc_matrix : log fold change in matrix form
#' @param entity scientific name of organism
#'
#' @import gage
#' @import utils
#' @import gage
#' @import pathview
#' @return nothing returned
#'

run_gage <- function(gsets, work.dir, same.dir,compare ,cpd_data, entity,
                     gene_data ,cpd.idtype = "kegg", gene_id_type ="entrez") {
    qcut<-0.1
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
    #    gage.out=gage.out[order(pgs.gene),]
    gage.out<- gage.out[order(pgs.gene, -gage.out[,"set.size"]),]
    write.table(gage.out, file = file.path(work.dir , paste0("gage.out.", 
    datatype , ".tsv")), sep="\t", col.names=NA, quote = FALSE)
    sig.i<- gage.out[,"q.val"]<qcut & !is.na(gage.out[,"q.val"])
    nsig<-sum(sig.i, na.rm=T)
    if(nsig>0) {
        gage.out.sig<-data.frame(gage.out)[sig.i,]
        ord1<-order(gage.out.sig[,"stat.mean"], decreasing=TRUE)
        gage.out.sig<-gage.out.sig[ord1,]
        write.table(gage.out.sig, file =file.path(work.dir, 
        paste0("gage.out.sig.", datatype , ".tsv")), sep="\t", 
                    col.names=NA, quote = FALSE)
        gpath_ids=rownames(gage.out.sig)
    } else {
        print("No gene set selected in 1d-test, select the top 3 instead!")
        #    grg=gage.res$greater
        #    grl=gage.res$less
        #    grg=grg[order(grg[,"p.val"], -grg[,"set.size"]),]
        #    grl=grl[order(grl[,"p.val"], -grl[,"set.size"]),]
        #    gpath_ids=unique(c(rownames(grg)[1:3],rownames(grl)[1:3]))
        gsel<-gage.out[,"set.size"]>0
        if(sum(gsel)>0)  gpath_ids<-rownames(gage.out)[1:min(sum(gsel),3)]
         }# gpath_ids=rownames(gage.out)[1:3]
    gpath_ids <- substr(gpath_ids, start = 4, stop = 8)
    print(gpath_ids)
    if (same.dir == FALSE){
        plotpathways(gage.dir = work.dir, entity, gpath_ids[1:6] , 
                     cpd_data = cpd_data, gene_data =gene_data ,
                     cpd.idtype = "kegg", gene_id_type ="entrez")}
        return(list("pathways_selected"= gpath_ids,"pgs.gene" = pgs.gene,
                "gage_result"= gage.out,"data_used" =fc_matrix , 
                "gene_sets"=gsets ))
}



#####
#gage.dir,entity,path.ids, fc_matrix,cpd_data = NULL