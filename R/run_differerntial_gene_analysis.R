#' Prepare to run standard DESeq2 or edgeR for differential gene expression analysis and plot
#' volcano plots
#'
#'
#'
#' @param cnts :counts of gene as data.frame or filename of count of gene,one count per sample per gene , sample in row
#' @param outdir : directory to store results of deseq2
#' @param entity : scientific name of organism to convert gene to symbol
#' @param gcompare : if the experiment is paired or unpaired
#' @param npca number of genes to use for pca
#' @param nheatmap number of genes for heatmap
#' @param phenofile_res returned object from process_phenofile
#' @param diff.tool either DeSeq2 or edgeR
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

run_differerntial_gene_analysis <- function(cnts, phenofile_res, outdir, entity, gcompare, npca, nheatmap, diff.tool = "DeSeq2"){
    
    if(!inherits(cnts, "data.frame")){
        if (file.exists(cnts)){
            cnts <- read.table(cnts,header = TRUE, sep = "\t" )
        } else{
            message("make sure cnts is of class data.frame or path to existing gene count file")
        }
    }

    cnts[cnts == 0] <- NA
    cnts <- as.data.frame(na.omit(cnts))    
    
    #find the two group
    #phenofile_res <- process_phenofile(phenofile)
    coldata <- phenofile_res$coldata
    SampleName <- phenofile_res$SampleName
    #filenames <- phenofile_res$FileName
    paired_info <- as.factor(as.character(phenofile_res$paired_info))
    
    cnts <- cnts[, coldata$SampleName]
    if (all(coldata$SampleName == colnames(cnts))) { # if this then proceed
        ref <- which(coldata$Class == levels(as.factor(coldata$Class))[1])
        samp <- which(coldata$Class == levels(as.factor(coldata$Class))[2])
        grp.idx <- NULL
        grp.idx[ref] <- "reference"
        grp.idx[samp] <- "sample"
    } else {
        message("make sure pheno file have only samples analysed")
    }
    if (gcompare=="paired"){
        newSampleName <- as.factor(as.character(paired_info)) #make sure this is in accordance
        grp.idx <- as.factor(grp.idx)
        coldat <-  DataFrame(cbind(grp.idx, newSampleName))
        formula_string <- "~newSampleName+grp.idx"
    } else { 
        coldat <- DataFrame(grp.idx = factor(grp.idx))
        formula_string <- "~grp.idx" }
    formula_object <- as.formula(formula_string)



    if (diff.tool == "DESeq2"){
        gene_data <- run_deseq2(cnts, 
                        outdir, entity,  formula_object,coldat, npca=npca, nheatmap=nheatmap ) 
    } else {
        gene_data <- run_edgeR(cnts,  outdir,grp.idx,entity  )
    }

return(gene_data)
}
