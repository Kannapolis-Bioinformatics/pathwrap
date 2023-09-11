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
#' @importFrom graphics plot.new
#' @import EnhancedVolcano EnhancedVolcano
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @import DESeq2
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom ComplexHeatmap pheatmap
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
  names(deseq2.fc) <- rownames(deseq2_res)
  exp.fc <- deseq2.fc
  table(is.na(deseq2_res$padj))
  write.table(
    deseq2_res,
    file.path(deseq2.dir, "DESEQ2_logfoldchange.txt"),
    sep = "\t",
    col.names = NA,
    row.names = TRUE,
    quote = FALSE
  )
  plot.new()
  tiff(
    file.path(deseq2.dir, "Volcano_deseq2.tiff"),
    units = "in",
    width = 15,
    height = 15,
    res = 300
  )
  plot(
    EnhancedVolcano::EnhancedVolcano(
      deseq2_res,
      x = "log2FoldChange",
      y = "pvalue",
      lab = rownames(deseq2_res)
    )
  )
  dev.off()
  #####################################
  sigs <- na.omit(deseq2_res)
  df <- as.data.frame(sigs)

  df.top <- df[(df$padj < 0.05) & (abs(df$log2FoldChange) > 2), ]
  if (dim(df.top)[1] > 19) {
    df.top <- na.omit(df.top[order(df.top$log2FoldChange,
      decreasing = TRUE )[seq_len(20)], ])
    ######################################
    message("Principle Componenet Analysis using VST from DESeq2")
    a <- try({
      vsd <- vst(dds, blind = TRUE, nsub =dim(df.top)[1] )
    })
    loadError <- (is(a, "try-error") | is(a, "error"))
    if (loadError == T) {
      message("transformation using varianceStabilizingTransformation")
      vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
    }
    mat <- assay(vsd)[rownames(df.top), colnames(cnts)]

    message("Now we are plotting PCA")
    plot.new()
    tiff(
      file.path(aligned_bam, "PCA_vst.tiff"),
      units = "in",
      width = 15,
      height = 15,
      res = 300
    )
    g <- plotPCA(vsd, intgroup = c("grp"))
    plot(g)
    dev.off()
    plot.new()

    message(
      "Also plotting heatmap of vst count of top 20 genes with
                    LFC more than 2 and padj less than 0.05"
    )
    tiff(
      file.path(aligned_bam, "heatmap_vst.tiff"),
      units = "in",
      width = 15,
      height = 15,
      res = 300
    )
    g <- pheatmap(as.numeric(mat), scale = "row", 
                 cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)
    plot(g)
    dev.off()
  }
  return(exp.fc)
}
