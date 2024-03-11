#' Wrapper for RNASeq data analysis from raw reads to pathway visualization
#'
#' wrapper that does quality control analysis of raw files, performs adapter and
#' quality trimming, builds genome index and does alignment , counts genes
#' perfoms differential gene expression analysis,does gene enrichment test using
#' GAGE and visualize the enriched pathways using pathview all using one
#' wrapper function. It has the ability to continue the analysis if it is
#' halted at any stage and generate quality pictures and generate comprehensive
#' analysis of the data.
#'
#' @param ref.dir : path to reference directory which contain
#' reference file(*.fa) and annotation file(*.gtf), can be NULL
#' @param phenofile : path to phenofile ; see note on test data to
#' see the format of phenofile
#' @param outdir : give the name of parent result dir , this can be existing or
#' not , rest of directory will be formed by program for organization
#' @param entity : Scientific name of species whose RNA is being analyzed
#' @param corenum : number of cores available for analysis #defaut 2
#' @param diff.tool : what differential tool to to use,
#' “DESEQ2” or “edgeR” available
#' @param compare : what is sample/experimental design you have,
#'    paired or unpaired, as.group
#' @param aligner : One of "Rhisat2" or "Rbowtie2"; Rbowtie2 can be very slow
#'    for human and eukaryotic species
#' @param keep_tmp : set TRUE if keeping the aligned bam files, if set FALSE,
#'    bam files are deleted
#' @param startover : if TRUE, all previous analysis is deleted 
#' if TRUE analysis starts from first step
#' @param cacheDir : directory where temporary files created during alignment
#' @param gene_id : gene id type 
#' @param cpd_id : compound id type
#' @param csamp  : sample column of compound
#' @param cref  : reference column of compound
#' @param mode :how to select pathways,"auto","combined"=using both cpd and gene
#' @param pid  kegg pathway ids
#' @param cdatapath path for the tsv file of compound
#' @param ccompare what is the comparision for sample and reference
#' @return the analysis status message
#'
#' @importFrom stringr str_replace_all
#' @import parallel
#' @import utils
#'
#' @example inst/example.R
#'
#' @export
pathviewwrap <- function(ref.dir = NA, phenofile = NA, outdir = "results",
                        entity = "Mus musculus",
                        corenum = 2, compare = "unpaired",
                        diff.tool = "DESeq2",keep_tmp = FALSE, 
                        startover = FALSE, cacheDir = NULL,
                        aligner = "Rhisat2", gene_id = NULL, 
                        cpd_id = "KEGG COMPOUND accession",csamp= NULL, cref = NULL,
                        mode="auto",ccompare=NA, pid= NULL, cdatapath=NA ) {
    on.exit(closeAllConnections())

    dirlist <- createdir(pos =1, outdir, entity, startover, keep_tmp)
    if (is.null(dirlist)){
        return(paste0("Please make sure the outdir is empty or start ",
                        "analysis with startover= FALSE", collapse="")) 
    }
    
    reference_paths <- sanity_check(ref.dir, outdir, entity, corenum, compare )
    
    if (is.null(reference_paths)) {
        if (is.na(ref.dir)){
            message("Please install the reference package")
        } else{
        message("Please make sure the file is bgzipped or ref.dir is writtable")
        }
        return("Please startover analysis with startover = TRUE")
    }

    genomeFile <- reference_paths[1]
    message("this is genome used")
    message(genomeFile)
    geneAnnotation <- reference_paths[2]
    message("this is annotation used")
    message(geneAnnotation)
    qc.dir <- dirlist[1]
    trim.dir <- dirlist[2]
    deseq2.dir <- dirlist[3]
    edger.dir <- dirlist[4]
    gage.dir <- dirlist[5]
    cset_dir <- dirlist[6]
    combined_dir <- dirlist[7]
    
    if (!file.exists(phenofile)) { ### TO DO make sure reference is first ANum
        message("Please provide phenofile with Class information")
    }
    coldata <- read.table(phenofile, sep = "\t", header = TRUE)
    if (colnames(coldata)[ncol(coldata)] != "Class") {
        message("Please make sure class information is in last column with
                colname 'Class'. ")
    }
    coldata$Class <- as.factor(coldata$Class)
    SampleName <- coldata$Sample
    filenames <- as.data.frame(coldata[, -c(1, ncol(coldata))])
    
    if (dim(filenames)[2] == 1) {
        endness <- "SE"
        fq.dir <- dirname(filenames[1, 1])
    } else if (dim(filenames)[2] == 2) {
        endness <- "PE"
        fq.dir <- dirname(filenames$FileName1[1])
    }
    
    if (!file.exists(file.path(qc.dir, "qc_heatmap.tiff",
                            fsep = .Platform$file.sep))) {
        message("STEP 1 ; running fastqc")
        message("this is qc.dir")
        message(qc.dir)
        res.fastqc <- run_qc(fq.dir, qc.dir, corenum)
        if (is.null(res.fastqc)) {
            return("Please make sure fastqc is available in system and
                accessible to R")
        }
    }
    if (length(list.files(trim.dir, "html"))<length(SampleName)){
    
        # just in case there is random component in run_fastp
        RNGkind("L'Ecuyer-CMRG")
        for (idxval in seq_len(length(SampleName))){
            run_fastp(SampleName[idxval], filenames[idxval,], 
                    endness, trim.dir, corenum)
        }
    }
    sampleFile <- writesampleFile(outdir, filenames,
            SampleName, trim.dir, endness)
    # make txdb from annotation
    
    txdbfilename <- paste0(gsub(" ", "", entity), "_txdbobj", collapse = "")
    if (!file.exists(file.path(outdir, txdbfilename,
                                fsep = .Platform$file.sep))) {
        message("STEP 2; making txdb obj")
        txdb <- make_txdbobj(
            geneAnnotation, corenum,
            genomeFile, entity, outdir
        )
    } else {
        txdb <- AnnotationDbi::loadDb(file.path(outdir, txdbfilename,
                                                fsep = .Platform$file.sep))
    }
    if (!file.exists(file.path(aligned_bam, "alltrimmedalignedobj.RDS",
                            fsep = .Platform$file.sep))) {
        message("STEP 3 : aligning the sequence")
        aligned_proj <- run_qAlign(
            corenum, endness, sampleFile, genomeFile,
            geneAnnotation, ref.dir, cacheDir, aligner
        )
    } else {
        aligned_proj <- readRDS(file.path(
            aligned_bam, "alltrimmedalignedobj.RDS",
            fsep = .Platform$file.sep  ))
    }
    if (!file.exists(file.path(outdir, "combinedcount.trimmed.RDS",
                            fsep = .Platform$file.sep))) {
        message("STEP 4: counting aligned sequences")
        cnts <- run_qCount(aligned_proj, corenum, outdir, txdb, entity)
    } else {
        cnts <- as.data.frame(readRDS(file.path(
            outdir, "combinedcount.trimmed.RDS", fsep = .Platform$file.sep
        )))
    }
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
    if (keep_tmp == FALSE) {
        message("deleting aligned bam files, bam file index and log files")
        unlink(list.files(file.path(outdir, "aligned_bam", 
        fsep = .Platform$file.sep), pattern = ".bam$|.bai$", 
        full.names = TRUE
        ))
    }
    if (!file.exists(file.path(deseq2.dir, "Volcano_deseq2.tiff",
                            fsep = .Platform$file.sep))) {
        message("STEP 5a ; running differential analysis using DESeq2")
        exp.fcncnts.deseq2 <- run_deseq2(cnts, grp.idx, deseq2.dir, entity)
    } else {
            deseq2.res.df <- read.table(
            file.path(deseq2.dir, "DESEQ2_logfoldchange.txt", 
                    fsep = .Platform$file.sep ),
            header = TRUE, sep = "\t", row.names = 1 )
        # works with gage
        exp.fcncnts.deseq2 <- deseq2.res.df$log2FoldChange
        names(exp.fcncnts.deseq2) <- rownames(deseq2.res.df)
    }
    if (!file.exists(file.path(edger.dir, "Volcano_edgeR.tiff",
                            fsep = .Platform$file.sep))) {
        message("STEP 5b ; running differential analysis using edgeR")
        exp.fcncnts.edger <- run_edgeR(cnts, grp.idx, edger.dir)
    } else {
        edger.res.df <- read.table(
            file.path( edger.dir, "edgeR_logfoldchange.txt",
                    fsep = .Platform$file.sep),
            header = TRUE, sep = "\t", row.names = 1 )
        # works with gage
        exp.fcncnts.edger <- edger.res.df$logFC
        names(exp.fcncnts.edger) <- rownames(edger.res.df)
    }
    if (diff.tool == "DESeq2") {
        exp.fc <- exp.fcncnts.deseq2
    } else {
        exp.fc <- exp.fcncnts.edger
    }
    if (!file.exists(file.path(gage.dir , "KEGG.sig.txt",
                            fsep = .Platform$file.sep))) {
        message("STEP 7 : running gene pathway analysis using GAGE")
        message(paste0(compare, "this is from pathviewwrap", collapse = ""))
        res_gage_gene <-run_pathway(entity, exp.fc, compare, gage.dir, 
                                    cnts, grp.idx)
        #return transfer properly TO DO
        gpath_ids <- res_gage_gene[1]
        pgs.gene <- res_gage_gene[2]
        gage_out <- res_gage_gene[3]
        gsets <- res_gage_gene[4]
    }
    if(!file.exists(file.path(cset_dir, "CKEGG.sig.txt",
                            fsep = .Platform$file.sep )) & !is.na(cdatapath)){
        message("STEP 8: running compound set analysis using GAGE")
        res_gage_cpd <-run_cpathway(cdatapath,cpd_id, csamp,cref, 
                                    ccompare,cset_dir )
        cpath_ids <- res_gage_cpd[1]
        pgs_cpd <- res_gage_cpd[2]
        gage_out_cpd <- res_gage_cpd[3]
        cpd_data <- res_gage_cpd[4]
    }
    if (mode == "combined"){
        message("STEP 9: running combined gene set analysis using GAGE")
        qcut <- 0.2
        path.ids <- run_combinedpath_analysis(gpath_ids, cpath_ids,gsets, 
                    pgs.gene,pgs_cpd, cset_dir, gage_out, gage_out_cpd, qcut)
        plotpathways(combined_dir,entity,path.ids, 
                    exp.fc,cpd_data = cpd_data)
    }
    onexistcleanup(ref.dir, entity)
    return("The analysis is complete")
}
