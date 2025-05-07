#' Run standard edgeR2 for differential gene expression analysis and
#' plot volcano plots
#'
#' run_edgeR takes counts and the list indicating reference and samples and
#' the directory where the results are stored and performs the edgeR analysis
#' The output is result table with columns of genes and logFC from result of
#' edgeR analysis and a volcanoplot.
#' It returns the logFC to be used by GAGE for gene set analysis.
#'
#'
#' @param cnts : counts of genes
#' @param grp.idx : index of the reference and sample for differential analysis as factor
#' @param outdir the  parent directory in which edgeR results will be stored
#' @import edgeR
#'
#' @return log fold expression values
#'

run_edgeR <- function(cnts, outdir, grp.idx){ 
    edger.dir <- file.path(outdir, "differential_analysis", "edgeR")
    if(!inherits(cnts)=="data.frame"){
        if (file.exists(cnts)){
            cnts <- read.table(cnts,header = TRUE, sep = "\t" )
        } else{
            message("make sure cnts in of class data.frame or path to existing gene count file")
        }
    }
    cnts[cnts == 0] <- NA
    cnts <- as.data.frame(na.omit(cnts))
    
    #volcano plot file
    edgeRvolcanoplot_filepath<- file.path(edger.dir, "Volcano_edgeR.tiff",fsep = .Platform$file.sep)
    if (!file.exists(edgeRvolcanoplot_filepath)){
        #do the edger
    
        dgel <- edgeR::DGEList(counts = cnts, group = factor(grp.idx))
        dgel <- edgeR::calcNormFactors(dgel)
        dgel <- edgeR::estimateCommonDisp(dgel)
        dgel <- edgeR::estimateTagwiseDisp(dgel)
        et <- edgeR::exactTest(dgel)
        edger.fc <- et$table$logFC
        names(edger.fc) <- rownames(et$table)
        exp.fc <- edger.fc
        write.table(
            et,
            file.path(edger.dir, "edgeR_logfoldchange.txt",
                    fsep = .Platform$file.sep),
            sep = "\t",
            col.names = TRUE,
            row.names = TRUE,
            quote = FALSE
        )
        tiff(
            edgeRvolcanoplot_filepath,
            units = "in",
            width = 15,
            height = 15,
            res = 300
        )
        plot(
            EnhancedVolcano::EnhancedVolcano(
                et$table,
                x = "logFC",
                y = "PValue",
                lab = rownames(et$table)
            )
        )
        dev.off()
    }
    return(list("logfoldchange"=exp.fc))

    ######
}
