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
#' @param entity : scientific name of organism to convert gene to symbol
#' @param SampleName list of the sample name to be used in paired experiment
#' @param compare : if the experiment is paired or unpaired
#' @import EnhancedVolcano EnhancedVolcano
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @import DESeq2
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom ComplexHeatmap pheatmap
#' @import pathview
#' @import gage
#' @return fold change values
#'

run_deseq2 <- function(cnts, grp.idx, deseq2.dir, entity,SampleName=NULL,
                    compare) {
    
    if (compare=="paired"){
        newSampleName <- as.factor(SampleName[[1]])
        grp.idx <- as.factor(grp.idx)
        coldat <-  DataFrame(cbind(grp.idx, newSampleName))
        formula_string <- "~newSampleName+grp.idx"
    } else { 
        coldat <- DataFrame(grp = factor(grp.idx))
        formula_string <- "~grp" }
    formula_object <- as.formula(formula_string)
    dds <- DESeqDataSetFromMatrix(cnts, colData = coldat, 
                            design = formula_object )
    dds <- DESeq(dds)
    deseq2_res <- results(dds)
    # direction of fc, depends on levels(coldat$grp), the first level
    # taken as reference (or control) and the second one as experiment.
    deseq2.fc <- deseq2_res$log2FoldChange
    names(deseq2.fc) <- rownames(deseq2_res)
    exp.fc <- deseq2.fc
    deseq2_res_ord <- deseq2_res[order(deseq2_res$pvalue),]
    write.table(
        deseq2_res_ord,
        file.path(deseq2.dir, "DESEQ2_logfoldchange.txt",
                fsep = .Platform$file.sep),
        sep = "\t", col.names = NA,     row.names = TRUE,    quote = FALSE)   
    data(bods, package = "gage", envir = environment())
    data(korg, package = "pathview", envir = environment())
    data(gene.idtype.list, package = "pathview", envir = environment())
    genesymbols <- eg2id(eg=rownames(deseq2_res), 
                        category = c("SYMBOL", "GENENAME"),
        org = unname(bods[bods[,3]==unname(korg[korg[,4]==entity, 3]) , 2]))
    tiff(
        file.path(deseq2.dir, "Volcano_deseq2.tiff",fsep = .Platform$file.sep),
        units = "in", width = 15,height = 15, res = 300)
    #plot has ensembl/gencode geneids
    plot(
        EnhancedVolcano::EnhancedVolcano(deseq2_res,
            x = "log2FoldChange", y = "pvalue",lab = genesymbols[,2]))
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
        if (dim(df.top)[1] > 200) {
            message("Now we are plotting PCA")
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
