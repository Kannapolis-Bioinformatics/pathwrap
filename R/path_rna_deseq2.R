#' Run standard DESeq2 for differential gene expression analysis and plot
#' volcano plots
#'
#' run_deseq2 takes counts and the list indicating reference and samples and the
#' directory where the results are stored and performs the deseq2 analysis
#' The output is result table with columns of genes and log2FoldChange from
#' result of deseq2 analysis and a volcanoplot.
#' For a paired experiment, sample name of the pair should be same but their class is different
#' It returns the log2foldChange to be used by GAGE for gene set analysis.
#'
#' @param cnts :counts of gene as data.frame or filename of count of gene,one count per sample per gene , sample in row
#' @param outdir : directory to store results of deseq2
#' @param entity : scientific name of organism to convert gene to symbol
#' @param coldat the dataframe for phenotype information of all sample
#' @param npca number of genes to use for pca
#' @param nheatmap number of genes for heatmap
#' @param formula_object string denoting formula
#' @import EnhancedVolcano EnhancedVolcano
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @import DESeq2
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom ComplexHeatmap pheatmap
#' @import pathview
#' @import gage
#' @return fold change values
#' @export 

run_deseq2 <- function(cnts, outdir, entity,  formula_object,coldat,npca, nheatmap){
    
    
    deseq2.dir <- file.path(outdir, "differential_analysis" , "DESeq2")
    kegg.gs.species <- kegg.gsets(entity)
    orgcode <- kegg.species.code(entity)
    volcanoplotfilepath <- file.path(deseq2.dir, "Volcano_deseq2.tiff",fsep = .Platform$file.sep)
    if(!file.exists(volcanoplotfilepath)){
        dds <- DESeqDataSetFromMatrix(cnts, colData = coldat, 
                                      design = formula_object )
        dds <- DESeq(dds)
        norm_counts<-counts(dds,normalized=TRUE)
        write.table(norm_counts, file= file.path(deseq2.dir, "normalized_count.txt"),col.names = T, row.names = T, quote=FALSE)
        deseq2_res <- results(dds)
        deseq2_res_copy <- deseq2_res
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
        
        
        # if(!all(rownames(cnts)%in% unlist(unname(kegg.gs.species$kg.sets))))
        # { #check if the use of "all" is appropriate
        
        # converting to entrez # what if gene id is not ensembl and what if
        # arabidopsis thaliana id.map might be ath or else thing
        rownames(deseq2_res_ord) <- str_remove(rownames(deseq2_res_ord),pattern= "\\.\\d+")
        rownames(deseq2_res) <- str_remove(rownames(deseq2_res),pattern= "\\.\\d+")
        
        #gene symbol in volcano plot
        #logfoldchange should have entrez id
        data(bods, package = "gage", envir = environment())
        data(korg, package = "pathview", envir = environment())
        data(gene.idtype.list, package = "pathview", envir = environment())
        rownames(deseq2_res) <- str_remove(rownames(deseq2_res), "\\.[0-9]+$")
        if (sum(rownames(deseq2_res) %in% unlist(unname(kegg.gs.species$kg.sets))) < 10) {
            org<-  unname(korg[korg[,4]==entity, 3])
            twoletter = paste0(toupper(substr(org, 1, 1)), substr(org, 2, 2))
            entrezid <- id2eg( ids = rownames(deseq2_res),category="ENSEMBL",org=twoletter ) 
            
            
        #kegg ~entrez
            genesymbols <- eg2id(eg=entrezid[,2], category = c( "ENTREZID","SYMBOL"),
                             org = unname(bods[bods[,3]==unname(korg[korg[,4]==entity, 3]) , 2]), unique.map=TRUE, na.rm=TRUE, keep.order=TRUE)
        
            for (val in 1:nrow(entrezid)){
                if(!is.na(entrezid[val,2])){
                    rownames(deseq2_res)[val] <- entrezid[val,2]
                }
            }
            labstoplot<- genesymbols[,2]
        }
        else{
            genesymbols <- eg2id(eg=rownames(deseq2_res), category = c( "ENTREZID","SYMBOL"),
                                 org = unname(bods[bods[,3]==unname(korg[korg[,4]==entity, 3]) , 2]), unique.map=TRUE, na.rm=TRUE, keep.order=TRUE)
            labstoplot<- genesymbols[,2]
        }
        deseq2.fc <- deseq2_res$log2FoldChange
        names(deseq2.fc) <- rownames(deseq2_res)
        exp.fc <- deseq2.fc
        write.table(deseq2_res, file.path(deseq2.dir, "DESEQ2_logfoldchange.txt",
                                          fsep = .Platform$file.sep), row.names=T, col.names  =T, sep ="\t")
        
        tiff(volcanoplotfilepath,
             units = "in", width = 15,height = 15, res = 300)
        #plot has ensembl/gencode geneids
        plot(
            EnhancedVolcano::EnhancedVolcano(deseq2_res,
                                             x = "log2FoldChange", y = "pvalue",lab = labstoplot))
        dev.off()
        
        
        vsd_matrix <- plotdeseqheatmap(deseq2_res_copy,dds,deseq2.dir, npca, nheatmap)
        message("plot return function 2 after retun")
        #write.table(vsd_matrix, file= file.path(deseq2.dir, "vsd_count.txt"),col.names = T, row.names = T, quote=FALSE)
        
    } else {
        exp_table <- read.table( file.path(deseq2.dir, "DESEQ2_logfoldchange.txt",
                                        fsep = .Platform$file.sep), row.names=1, header =T, sep = "\t")
       exp.fc <-  exp_table[,2]
       names(exp.fc) <- rownames(exp_table)
       norm_counts <- data.matrix(read.table(file.path(deseq2.dir, "normalized_count.txt"),header = T, row.names = 1, quote=""))
       print("this isthe size of log fold change being read")
       print(dim(exp.fc))
       #vsd_matrix <- read.table(file.path(deseq2.dir, "vsd_count.txt"),header = T, row.names = 1, quote="")
    }
    
    return(list("logfoldchange"=exp.fc, "normalized_count"=norm_counts, "keggorgcode" = orgcode))#, "vsd_count"= vsd_matrix))
}

#' plot the result of result of deseq2    
#' @param deseq2_res result of deseq2
#' @param dds deseq2 object
#' @param deseq2.dir directory where aligned bam are stored
#' @param npca number of genes to use for pca
#' @param nheatmap number of gens to be used for heatmap
#' @return vsd_matrix just returns    
plotdeseqheatmap <- function(deseq2_res,dds,deseq2.dir,npca, nheatmap){
    #####################################
    dds$grp.idx <- as.factor(dds$grp.idx)
    sigs <- na.omit(deseq2_res)
    df <- as.data.frame(sigs)
    df.top <- df[(df$padj < 0.05) & (abs(df$log2FoldChange) > 2), ]
    message("total number of differentially expressed gene, padj<0.05 and logfoldchange>2 is  ", dim(df.top)[1])
    vsd_matrix<-NULL
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
    vsd_matrix<- assay(vsd)
    
    if (dim(df.top)[1] >npca) {
        message("Now we are plotting PCA")
        tiff(
            file.path(deseq2.dir, "PCA_vst.tiff",fsep = .Platform$file.sep),
            units = "in",  width = dim(df.top)[2], height = dim(df.top)[2],  res = 300)
        
        g <- plotPCA(vsd, intgroup = c("grp.idx")) +
            theme(
                axis.text.x = element_text(size = 14),     # x-axis tick labels
                axis.title.x = element_text(size = 20),   # x-axis title
                axis.text.y  = element_text(size = 14),     # x-axis tick labels
                axis.title.y = element_text(size = 20)   # x-axis title
            )
        
        plot(g)
        dev.off()
    }
    df.top <- na.omit(df.top[order(df.top$log2FoldChange,
                                   decreasing = TRUE )[seq_len(20)], ])
    if (dim(df.top)[1] > nheatmap) {
        rowstoselect <- match(rownames(df.top)[seq_len(nheatmap+1)], 
                              rownames(assay(vsd)))
        mat <- assay(vsd)[rowstoselect , ]  
        message(
            "Also plotting heatmap of vst count of top 20 genes with
                LFC more than 2 and padj less than 0.05"
        )
        tiff(
            file.path(deseq2.dir, "heatmap_vst.tiff",fsep = .Platform$file.sep),
            units = "in",  width = 15, height = 15,res = 300)
        g <- pheatmap(as.matrix(mat), scale = "row", cluster_rows=FALSE, 
                      show_rownames=TRUE,cluster_cols=FALSE,fontsize = 15, 
                      heatmap_legend_param = list(column_names_rot=45))
            
        plot(g)#+ theme(axis.text.x = element_text(angle = 45, vjust = 1))
        dev.off()
    }
    message("plot return function vsd_matrix")
    return(vsd_matrix)
}
