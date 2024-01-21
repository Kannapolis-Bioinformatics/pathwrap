#' Run standard DESeq2 for differential gene expression analysis and plot
#' volcano plots
#'
#' run_deseq2 takes counts and the list indicating reference and samples and the
#' directory where the results are stored and performs the deseq2 analysis
#' The output is result table with columns of genes and log2FoldChange from
#' result of deseq2 analysis and a volcanoplot.
#' It returns the log2foldChange to be used by GAGE for gene set analysis.
#'
#' @param cnts :counts of gene
#' @param grp.idx : index of reference and samples for differential analysis
#' @param deseq2.dir : directory to store results of deseq2
#' @import EnhancedVolcano EnhancedVolcano
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @import DESeq2
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom ComplexHeatmap pheatmap
#' @import pathview
#' @return fold change values
#'

run_deseq2 <- function(cnts, grp.idx, deseq2.dir) {
    coldat <- DataFrame(grp = factor(grp.idx))
    dds <-
        DESeqDataSetFromMatrix(cnts, colData = coldat, design = ~grp)
    dds <- DESeq(dds)
    deseq2_res <- results(dds)
    # direction of fc, depends on levels(coldat$grp), the first level
    # taken as reference (or control) and the second one as experiment.
    deseq2.fc <- deseq2_res$log2FoldChange
    #names(deseq2.fc) <- rownames(deseq2_res)
    genesymbols <- eg2id(eg=rownames(deseq2_res), 
        org = unname(bods[bods[,3]==unname(korg[korg[,4]==entity, 3]) , 2])
    names(deseq2.fc) <-genesymbols 
    exp.fc <- deseq2.fc
    table(is.na(deseq2_res$padj))
    #org = common_name of organism from bods
    genesymbols <- eg2id(eg=rownames(deseq2_res), 
    org = unname(bods[bods[,3]==unname(korg[korg[,4]==entity, 3]) , 2])
    write.table(
        deseq2_res,
        file.path(deseq2.dir, "DESEQ2_logfoldchange.txt",
                fsep = .Platform$file.sep),
        sep = "\t", col.names = NA,     row.names = TRUE,    quote = FALSE)
    tiff(
        file.path(deseq2.dir, "Volcano_deseq2.tiff",fsep = .Platform$file.sep),
        units = "in", width = 15,height = 15, res = 300)
    #plot has ensembl/gencode geneids
    plot(
        EnhancedVolcano::EnhancedVolcano(deseq2_res,
            x = "log2FoldChange", y = "pvalue",lab = rownames(deseq2_res)))
    dev.off()
    plotdeseqheatmap(deseq2_res,dds,deseq2.dir)
    return(exp.fc)
}

#' plot the result of result of deseq2    
#' @param deseq2_res result of deseq2
#' @param dds deseq2 object
#' @param deseq2.dir directory where aligned bam are stored
#' @return just returns    
plotdeseqheatmap <- function(deseq2_res,dds,deseq2.dir){
    #####################################
    sigs <- na.omit(deseq2_res)
    df <- as.data.frame(sigs)
    df.top <- df[(df$padj < 0.05) & (abs(df$log2FoldChange) > 2), ]
        ######################################
        message("Principle Componenet Analysis using VST from DESeq2")
        a <- try({
            vsd <- vst(dds, blind = TRUE, nsub =dim(df.top)[1] )
        }, silent = TRUE)
        loadError <- (is(a, "try-error") | is(a, "error"))
        if (loadError == TRUE) {
            message("transformation using varianceStabilizingTransformation")
            vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
        }
        message("Now we are plotting PCA")
        if (dim(df.top)[1] > 200) {
        tiff(
            file.path(deseq2.dir, "PCA_vst.tiff",fsep = .Platform$file.sep),
            units = "in",  width = 15, height = 15,  res = 300)
        g <- plotPCA(vsd, intgroup = c("grp"))
        plot(g)
        dev.off()
    }
        df.top <- na.omit(df.top[order(df.top$log2FoldChange,
                                    decreasing = TRUE )[seq_len(20)], ])
        if (dim(df.top)[1] > 19) {
        rowstoselect <- match(rownames(df.top)[seq_len(20)], 
                            rownames(assay(vsd)))
        mat <- assay(vsd)[rowstoselect , ]  
        message(
            "Also plotting heatmap of vst count of top 20 genes with
                                        LFC more than 2 and padj less than 0.05"
        )
        tiff(
            file.path(deseq2.dir, "heatmap_vst.tiff",fsep = .Platform$file.sep),
            units = "in",  width = 15, height = 15,res = 300)
        g <- pheatmap(as.matrix(mat), scale = "row", 
                    cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE)
        plot(g)
        dev.off()
    }
    return(invisible(NULL))
}
