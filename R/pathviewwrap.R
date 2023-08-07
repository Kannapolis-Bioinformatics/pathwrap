#'  Wrapper for RNASeq data analysis from raw reads to pathway visualization
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
#' @param outdir  : give the name of parent result dir , this can be existing or
#' not , rest of directory will be formed by program for organization
#' @param entity  : Scientific name of species whose RNA is being analyzed
#' @param corenum : number of cores available for analysis #defaut 2
#' @param diff.tool :  what differential tool to to use,
#' “DESEQ2” or “edgeR” available
#' @param compare : what is sample/experimental design you have,
#'  paired or unpaired, as.group
#' @param seq_tech : Illumina, pacbio or nanopore
#' @param aligner : One of "Rhisat2" or "Rbowtie2"; Rbowtie2 can be very slow
#'  for human and eukaryotic species
#' @param keep_tmp : set TRUE if keeping the aligned bam files, if set FALSE,
#'  bam files are deleted
#' @param rerun : if FALSE the previously complete step will not be rerunned,
#' if TRUE analysis starts from first step
#' @param cacheDir : directory where temporary files created during alignment
#' are stored
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
                        diff.tool = "DESeq2", seq_tech = "Illumina",
                        keep_tmp = FALSE, rerun = FALSE, cacheDir = NULL,
                        aligner) {
    on.exit(closeAllConnections())
    aligned_bam <- NA
    dirlist <- sanity_check(ref.dir,
        pos = 1, outdir, entity, corenum, compare,
        rerun
    )
    if (is.null(dirlist)) {
        message("Please install the reference package")
        return("Please rerun analysis with rerun = TRUE")
    }
    qc.dir <- dirlist[1]
    trim.dir <- dirlist[2]
    genomeFile <- dirlist[3]
    message("this is genome File")
    message(genomeFile)
    geneAnnotation <- dirlist[4]
    message("this is geneAnnotation")
    message(geneAnnotation)
    deseq2.dir <- dirlist[5]
    edger.dir <- dirlist[6]
    gage.dir <- dirlist[7]

    # duplicated codes
    # if (!file.exists(phenofile)){ ###TO DO make sure reference is first ANume
    #   message("Please provide phenofile with Class information")
    # }
    # coldata <- read.table(phenofile, sep = "\t", header = TRUE)
    # if(colnames(coldata)[ncol(coldata)]!="Class"){
    #   message("Please make sure class information is in last column with
    # colname 'Class' . ")
    # }
    # coldata$Class <- as.factor(coldata$Class)
    # SampleName <- coldata$Sample
    # filenames <- coldata[,-c(1,ncol(coldata))]

    # if(is.null(dim(filenames))){
    #   endness <- "SE"
    #   fq.dir <-  dirname(filenames[1])
    # } else if(dim(filenames)[2] == 2){
    #   endness <- "PE"
    #   fq.dir <- dirname(filenames$FileName1[1])
    #
    # }

    # run the fastqc


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

    if (!file.exists(file.path(qc.dir, "qc_heatmap.tiff"))) {
        message("STEP 1 ; running fastqc")
        message("this is qc.dir")
        message(qc.dir)
        res.fastqc <- run_qc(fq.dir, qc.dir, corenum)
        if (is.null(res.fastqc)) {
            return("Please make sure fastqc is available in system and
                   accessible to R")
        }
    }
    sampleFile <- file.path(outdir, "sampleFile.txt")
    rawfileName <- as.data.frame(vapply(filenames, function(x) basename(x), 
                character(dim(filenames)[1])))


    fastp_files_name <- as.data.frame(vapply(as.list(rawfileName),
            function(x) str_replace_all(x, ".fastq.gz$", "_trimmed.fastq.gz"),
            character(dim(filenames)[1])))
    FileName <- vapply(as.list(fastp_files_name),function(x) file.path(trim.dir, x),
                   character(dim(filenames)[1]) )
  
    if (endness == "SE") {
        write.table(
            file = sampleFile, sep = "\t",
            as.data.frame(cbind(FileName, SampleName)),
            col.names = c("FileName", "SampleName"),
            quote = FALSE, row.names = FALSE
        )
    } else {
        write.table(
            file = sampleFile, sep = "\t",
            as.data.frame(cbind(
                FileName[, 1],
                FileName[, 2], SampleName
            )),
            col.names = c("FileName1", "FileName2", "SampleName"),
            quote = FALSE, row.names = FALSE
        )
    }

    # just in case there is random component in run_fastp
    RNGkind("L'Ecuyer-CMRG")
    if (Sys.which("fastp") == "") {
        return("Please make sure fastp is available for R to use")
    }

    cl <- makeCluster(corenum)
    clusterExport(cl, c(
        "seq_tech", "endness", "FileName", "filenames",
        "trim.dir"
    ),
    envir = environment()
    ) # .GlobalEnv) ??
    parSapply(cl, SampleName, run_fastp)
    message("the trim run is complete")
    stopCluster(cl)

    # to check if all the nodes run fine
    # bad <- sapply(r, inherits, what = "try-error") # r<- mclappy()

    # make txdb from annotation
    txdbfilename <- paste0(gsub(" ", "", entity), "_txdbobj")
    if (!file.exists(file.path(outdir, txdbfilename))) {
        message("STEP 2; making txdb obj")
        txdb <- make_txdbobj(
            geneAnnotation, corenum,
            genomeFile, entity, outdir
        )
    } else {
        txdb <- AnnotationDbi::loadDb(paste0(
            outdir, "/", gsub(" ", "", entity),
            "_txdbobj"
        ))
    }

    if (!file.exists(file.path(aligned_bam, "alltrimmedalignedobj.RDS"))) {
        # setwd(outdir)
        message("STEP 3 : aligning the sequence")
        aligned_proj <- run_qAlign(
            corenum, endness, sampleFile, genomeFile,
            geneAnnotation, ref.dir, cacheDir, aligner
        )
    } else {
        aligned_proj <- readRDS(file.path(
            aligned_bam,
            "alltrimmedalignedobj.RDS"
        ))
    }

    if (!file.exists(file.path(outdir, "combinedcount.trimmed.RDS"))) {
        message("STEP 4: counting aligned sequences")
        cnts <- run_qCount(aligned_proj, corenum, outdir, txdb, entity)
    } else {
        cnts <- as.data.frame(readRDS(file.path(
            outdir, "combinedcount.trimmed.RDS"
        )))
    }
    # message("these are samplename for cnts , cnts[, coldata$SampleName] ")
    # message(coldata$SampleName)
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
        # unlink(file.path(outdir, "aligned_bam", "*bam*"))
        unlink(list.files(file.path(outdir, "aligned_bam"),
            pattern = ".bam$|.bai$", full.names = TRUE
        ))
    }
    deseq_volcano_plot <- paste0(deseq2.dir, "/Volcano_deseq2.tiff")
    if (!file.exists(deseq_volcano_plot)) {
        message("STEP 5a ; running differential analysis using DESeq2")
        exp.fcncnts.deseq2 <- run_deseq2(cnts, grp.idx, deseq2.dir)
    } else {
        deseq2.res.df <- read.table(
            file.path(
                deseq2.dir,
                "DESEQ2_logfoldchange.txt"
            ),
            header = TRUE, sep = "\t", row.names = 1
        )
        # works with gage
        exp.fcncnts.deseq2 <- deseq2.res.df$log2FoldChange
        names(exp.fcncnts.deseq2) <- rownames(deseq2.res.df)
    }
    edger_volcano_plot <- paste0(edger.dir, "Volcano_edgeR.tiff")
    if (!file.exists(edger_volcano_plot)) {
        message("STEP 5b ; running differential analysis using edgeR")
        exp.fcncnts.edger <- run_edgeR(cnts, grp.idx, edger.dir)
    } else {
        edger.res.df <- read.table(
            file.path(
                edger.dir, "edgeR_logfoldchange.txt"
            ),
            header = TRUE, sep = "\t", row.names = 1
        )
        # works with gage
        exp.fcncnts.deseq2 <- edger.res.df$log2FC
        names(exp.fcncnts.deseq2) <- rownames(edger.res.df)
    }

    # setwd(gage.dir)
    # chosing to use deseq2 result or edger result for gage
    if (diff.tool == "DESeq2") {
        exp.fc <- exp.fcncnts.deseq2
    } else {
        exp.fc <- exp.fcncnts.edger
    }
    if (!file.exists("*.txt")) {
        message("STEP 6 : running pathway analysis using GAGE")
        message(paste0(compare, "this is from pathviewwrap"))
        run_pathway(entity, exp.fc, compare, gage.dir, cnts, grp.idx)
    }
    return("The analysis is complete")
}

# when loading pathview??
##### Loading required namespace: org.Mm.eg.db
# Installing package(s) 'org.Mm.eg.db'
# installing the source package ‘org.Mm.eg.db’
