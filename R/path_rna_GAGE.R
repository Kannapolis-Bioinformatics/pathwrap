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
run_gage <- function(gsets, work.dir, same.dir,compare ,fc_matrix , entity) {
    qcut<-0.2
    gage.dir <- dirname(dirname(work.dir))
    fc.kegg.p <- gage( fc_matrix, gsets = gsets, ref = NULL,samp = NULL,
                    same.dir = same.dir, compare = compare)
    sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(
        fc.kegg.p$greater[, "q.val"])
    path.ids <- rownames(fc.kegg.p$greater)[sel]
    anla_type <- "KEGG"
    if (same.dir == TRUE) {
        anla_type <- "GO"
        gage.dir <- file.path(gage.dir , "GO",fsep = .Platform$file.sep)
        sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(
            fc.kegg.p$less[, "q.val"])
        path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
        write.table(fc.kegg.p$less,
                    sep = "\t", file = file.path(work.dir, paste0(
    "fc.", anla_type, ".p.less.txt"), fsep = .Platform$file.sep))
        path.ids <- c(path.ids[seq_len(3)], path.ids.l[seq_len(3)])
    }
    path.ids <- substr(path.ids, 1, 8)
    # write.table(fc.kegg.p$greater, sep = "\t",
    #             file = file.path(work.dir, paste0(
    # "fc.", anla_type,"p.greater.txt"), fsep = .Platform$file.sep))
    # visualize top 3 pathways
    if (same.dir == FALSE) { # run pathview only for KEGG pathways
        plotpathways(gage.dir,entity,path.ids, fc_matrix)
        }
    # kegg.sig <- sigGeneSet(fc.kegg.p,
    # outname = paste0(entity, anla_type, ".sig", basename(work.dir)
    #                     ), pdf.size = c(17, 17), heatmap = FALSE)
    # # wont give heatmap for fold change used in gage
    # write.table(kegg.sig$greater,
    #             file = file.path(gage.dir, paste0(anla_type, ".sig.txt"),
    #                             fsep = .Platform$file.sep),sep = "\t")
    gage.res <- fc.kegg.p 
    pgs=gage.res$greater[,"p.val"]
    pms=cbind(pgs*2, (1-pgs)*2)
    pgs.gene=apply(pms, 1, function(x) min(x))
    qgs.gene=p.adjust(pgs.gene, method = "BH")

    colnames(pms)=c("p.up", "p.dn")
    gage.out=cbind(gage.res$greater[, c(2,5)], pms/2, p.val=pgs.gene, q.val=qgs.gene)
    #    gage.out=gage.out[order(pgs.gene),]
    gage.out=gage.out[order(pgs.gene, -gage.out[,"set.size"]),]
    write.table(gage.out, file = file.path(gage.dir, 
                                           paste0(anla_type, ".res.tsv"),
     fsep = .Platform$file.sep), sep="\t", col.names=NA, quote = FALSE)



    ### significant.genesets
    sig.i=gage.out[,"q.val"]<qcut & !is.na(gage.out[,"q.val"])
    nsig=sum(sig.i, na.rm=T)
    if(nsig>0) {
        gage.out.sig=data.frame(gage.out)[sig.i,]
        ord1=order(gage.out.sig[,"stat.mean"], decreasing=T)
        gage.out.sig=gage.out.sig[ord1,]
        write.table(gage.out.sig, file =file.path(gage.dir, 
                    paste0(anla_type, ".res.sig.tsv"),fsep = .Platform$file.sep),sep="\t", 
                    col.names=NA, quote = FALSE)

        gpath.ids<-rownames(gage.out.sig)
        #pdf(file.path(gage.dir, 
                      #paste0(anla_type, ".gage.res.gs.heatmap.pdf"),fsep = .Platform$file.sep))
        #gage:::gs.heatmap(gage.res$stats[gpath.ids, -1], limit = 5, 
                          #main = "GAGE test statistics")
        #dev.off()

    } else {
        print("No gene set selected in 1d-test, select the top 3 instead!")
        #    grg=gage.res$greater
        #    grl=gage.res$less
        #    grg=grg[order(grg[,"p.val"], -grg[,"set.size"]),]
        #    grl=grl[order(grl[,"p.val"], -grl[,"set.size"]),]
        #    gpath.ids=unique(c(rownames(grg)[1:3],rownames(grl)[1:3]))
        gsel=gage.out[,"set.size"]>0
        if(sum(gsel)>0)  gpath.ids=rownames(gage.out)[1:min(sum(gsel),3)]
        # gpath.ids=rownames(gage.out)[1:3]
    }
    return(list("pathways"= gpath.ids,"pgs.gene" = pgs.gene,
                "gage_result"= gage.out,"data_used" =fc_matrix , 
                "gene_sets"=gsets ))
    
    #return(invisible(x=NULL))
}
