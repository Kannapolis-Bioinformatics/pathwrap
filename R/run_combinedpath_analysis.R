#' Run the gene and compound set analysis
#' runs gage for gene and compound data and combines the result using gamma function
#' 
#' @param outdir main directory for storing output of the process
#' @param gene_data logfold change data or cnts
#' @param entity Scientific name of species whose RNA/compound is being analyzed
#' @param gcompare how samples should be compared for gene data
#' @param gage.out.cpd_res the out object of cpound analysis
#' @param qcut threshold for pathway selection
#' @import gage
#' @return pathids important pathways of interest
#' @export
run_combinedpath_analysis <- function(outdir, gene_data, entity, gcompare,gage.out.cpd_res,qcut ){
#     work.dir<- file.path(outdir, "combined_analysis")
#     gsets<-  kegg.gsets()$kg.sets
#     gage.out_res <- run_gage2(gene_data=fc_matrix, gsets= kegg.gsets()$kg.sets, ref=NULL, samp=NULL, same.dir=TRUE, compare=compare, work.dir=work.dir,gene_id_type = "ENTREZ")
#     gpath_ids <- gage.out_res$pids 
#     cpath_ids <- gage.out.cpd_res$pids 
# 
#     
#     combinedpath_analysis(gpath_ids, cpath_ids,gsets, pgs.gene,
#                           pgs.cpd, somdir, gage.out, gage.out.cpd,qcut=0.01)
# }
    gsets <- kegg.gsets(entity)$kg.sets
    work.dir <- file.path(outdir, "gage_results", "combined_analysis_kegg")
    qcut<-0.01
    gage.dir <- dirname(dirname(work.dir))
    fc.kegg.p <- gage( gene_data, gsets = kegg.gsets(entity)$kg.sets, ref = NULL,samp = NULL,
                       same.dir = F, compare = gcompare)
    pgs<-fc.kegg.p$greater[,"p.val"]
    pms<-cbind(pgs*2, (1-pgs)*2)
    pgs.gene<-apply(pms, 1, function(x) min(x))
    qgs.gene<-p.adjust(pgs.gene, method = "BH")
    colnames(pms)<-c("p.up", "p.dn")
    gage.out<-cbind(fc.kegg.p$greater[, c(2,5)], pms/2, p.val=pgs.gene, 
                    q.val=qgs.gene)
    gage.out<- gage.out[order(pgs.gene, -gage.out[,"set.size"]),] #order based on min of p-value
    write.table(gage.out, file = file.path(work.dir , paste0("gage.out.", 
                                                             "gene" , ".tsv")), sep="\t", col.names=NA, quote = FALSE)
    sig.i<- gage.out[,"q.val"]<qcut & !is.na(gage.out[,"q.val"])
    nsig<-sum(sig.i, na.rm=TRUE)
    if(nsig>0) {
        gage.out.sig<-data.frame(gage.out)[sig.i,]
        ord1<-order(gage.out.sig[,"stat.mean"], decreasing=TRUE) #why?
        gage.out.sig<-gage.out.sig[ord1,]
        gpath_ids<-rownames(gage.out.sig)
    } else {
        message(paste0("No ", "gene",  " set selected in GAGE test, top 3 ", 
                       "gene",  "set plotted instead!", collapse = ""))
        gsel<-gage.out[,"set.size"]>0
        if(sum(gsel)>0) gpath_ids<-rownames(gage.out)[seq_len(min(sum(gsel),3))] }
    gpath_ids <- substr(gpath_ids, start = 4, stop = 8)
    
    
    cpath_ids<- gage.out.cpd_res$pids
    cgage.out <- gage.out.cpd_res$gage.out.cpd
    pgs.c<-cgage.out$greater[,"p.val"]
    pms.c<-cbind(pgs.c*2, (1-pgs.c)*2)
    pgs.cpd<-apply(pms.c, 1, function(x) min(x))
    qgs.cpd<-p.adjust(pgs.cpd, method = "BH")
    colnames(pms.c)<-c("p.up", "p.dn")
    gage.out.cpd<-cbind(cgage.out$greater[, c(2,5)], pms.c/2, p.val=pgs.cpd, 
                    q.val=qgs.cpd)
    
    
   pathids <- combinedpath_analysis(gpath_ids, cpath_ids,gsets, pgs.gene,
    pgs.cpd, work.dir, gage.out, gage.out.cpd,qcut=0.01)
   return(pathids)
}
   