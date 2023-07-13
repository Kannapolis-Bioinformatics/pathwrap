#' Title
#'
#' @param cnts : counts of genes
#' @param grp.idx : index of the reference and sample for differential analysis
#' @param edgeR.dir : the directory in which edgeR results will be stored
#'
#' @import edgeR
#'
#' @return exp.fc : log fold expression values
#' @export
#'

run_edgeR <- function(cnts,grp.idx, edger.dir){
    library(edgeR)
    dgel <- edgeR::DGEList(counts=cnts, group=factor(grp.idx))
    dgel <- edgeR::calcNormFactors(dgel)
    dgel <- edgeR::estimateCommonDisp(dgel)
    dgel <- edgeR::estimateTagwiseDisp(dgel)
    et <- edgeR::exactTest(dgel)
    edger.fc=et$table$logFC
    names(edger.fc)=rownames(et$table)
    exp.fc=edger.fc
    write.table(et , file.path(edger.dir, "edgeR_logfoldchange.txt"), sep = "\t", col.names =TRUE, row.names =TRUE, quote =FALSE)
    tiff(paste0(edger.dir, "/Volcano_edgeR.tiff"), units="in", width=15, height=15, res=300)
    plot(EnhancedVolcano::EnhancedVolcano(et$table, x ='logFC', y="PValue", lab=rownames(et$table)))
    dev.off()
    return(exp.fc)

  ######
}
