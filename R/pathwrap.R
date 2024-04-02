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
#' reference file(*.fa) and annotation file(*.gtf), can be NA
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
#' @param is.test to ask for prompt during test , is.test = TRUE, for test
#' @param ccompare what is the comparision for sample and reference
#' @return the analysis status message
#' @importFrom stringr str_replace_all
#' @import parallel
#' @import utils
#'
#' @example inst/example.R
#'
#' @export
pathwrap <- function(ref.dir = NA, phenofile = NA, outdir = "results",
                        entity = "Mus musculus",cref = NULL,
                        corenum = 2, compare = "unpaired",
                        diff.tool = "DESeq2",keep_tmp = FALSE, 
                        startover = FALSE, cacheDir = NULL,
                        aligner = "Rhisat2", gene_id = NULL, 
                        cpd_id = "KEGG",csamp= NULL, 
                        mode="auto",ccompare=NA, pid= NULL,
                        cdatapath=NA ,is.test=FALSE) {
    on.exit(closeAllConnections())
    setupobj<- setup_n_qc( outdir, entity, startover, keep_tmp,ref.dir, corenum,
                        compare, phenofile, is.test)
    filenames <- unname(setupobj["filenames"]) 
    coldata <- setupobj["phendata"][[1]]
    SampleName <- unname(setupobj["SampleName"] )
    endness <- unname(setupobj["endness"]) 
    genomeFile<- setupobj["reference_paths"][[1]][1]
    geneAnnotation<- setupobj["reference_paths"][[1]][2]
    trim.dir <- file.path(outdir, "fastp_results", fsep = .Platform$file.sep)
    prepobj<- prepare_for_alignment(outdir,filenames,SampleName,
                                    endness, entity, geneAnnotation)
    sampleFile <- file.path(outdir, "sampleFile.txt")
    aligned_bam <- file.path(outdir, "aligned_bam", fsep = .Platform$file.sep)
    alignobj <-align_n_count(aligned_bam,corenum, endness,keep_tmp , 
                            sampleFile,geneAnnotation, ref.dir=ref.dir,
                            cacheDir, aligner,outdir, prepobj, entity, 
                            compare, genomeFile,coldata)
    deseq2.dir <- file.path(outdir, "differential_analysis/DESeq2")
    
    diff_ana_obj<- diff_analysis(outdir,cnts=as.data.frame(alignobj$counts), 
                        grp.idx=alignobj$grp.idx[[1]], entity, 
                        SampleName,compare,diff.tool )
    genecptsetanalysis(outdir, entity,exp.fc=diff_ana_obj, compare, 
            cnts=as.data.frame(alignobj$counts), cdatapath,cpd_id, 
            csamp,cref,ccompare, mode, pid)
            onexistcleanup(ref.dir, entity)
    return("The analysis is complete")
}

#' This is an internal function of R which creates the result files and checks 
#' if all necessary conditions are meet to run the wrapper
#' @param ref.dir : path to reference directory which contain
#' reference file(*.fa) and annotation file(*.gtf), can be NULL

#' @param outdir : give the name of parent result dir , this can be existing or
#' not , rest of directory will be formed by program for organization
#' @param entity : Scientific name of species whose RNA is being analyzed
#' @param corenum : number of cores available for analysis #defaut 2
#' @param compare : what is sample/experimental design you have,
#'    paired or unpaired, as.group
#' @param keep_tmp : set TRUE if keeping the aligned bam files, if set FALSE,
#'    bam files are deleted
#' @param startover : if TRUE, all previous analysis is deleted 
#' if TRUE analysis starts from first step
#' @param phenofile the path to the phenotype file; format in vignettes
#' @param is.test to ask for prompt during test , is.test = TRUE, for test
#' @return returnobject with      
#' @noRd    
setup_n_qc  <- function( outdir, entity, startover, keep_tmp,ref.dir,
                        corenum, compare, phenofile, is.test){
    dirlist <- createdir(pos =1, outdir, entity, startover, keep_tmp,is.test)
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
    message("this is annotation used")
    if (!file.exists(phenofile)) { ### TO DO make sure reference is first ANum
        message("Please provide phenofile with Class information")}
    coldata <- read.table(phenofile, sep = "\t", header = TRUE)
    if (colnames(coldata)[ncol(coldata)] != "Class") {
        message("Please make sure class information is in last column with
                colname 'Class'. ")}
    coldata$Class <- as.factor(coldata$Class)
    SampleName <- coldata$SampleName
    filenames <- as.data.frame(coldata[, -c(1, ncol(coldata))])
    if (dim(filenames)[2] == 1) {
        endness <- "SE"
        fq.dir <- dirname(filenames[1, 1])
    } else if (dim(filenames)[2] == 2) {
        endness <- "PE"
        fq.dir <- dirname(filenames$FileName1[1]) }
    qc.dir <- file.path(outdir, "fastqc_results")
    if (!file.exists(file.path(qc.dir, "qc_heatmap.tiff",
                            fsep = .Platform$file.sep))) {
        message("STEP 1 ; running fastqc")
        res.fastqc <- run_qc(fq.dir, qc.dir, corenum)
        if (is.null(res.fastqc)) {
            return("Please make sure fastqc is available in system and
                accessible to R")
        } }
    if (length(list.files(dirlist[2], "html"))<length(SampleName)){
        RNGkind("L'Ecuyer-CMRG")
        for (idxval in seq_len(length(SampleName))){
            run_fastp(SampleName[idxval], filenames[idxval,], 
                    endness, dirlist[2], corenum)
        }   }
    return(list("filenames"= filenames, "SampleName" = SampleName, 
    "endness"= endness,"reference_paths" = reference_paths,"phendata"= coldata))
}

#' This function does prepares for the alignment
#' This function prepares for alignement
#' @param outdir : give the name of parent result dir , this can be existing or
#' not , rest of directory will be formed by program for organization
#' @param entity : Scientific name of species whose RNA is being analyzed
#' @param filenames dataframe of location of raw files
#' @param SampleName names of sample 
#' @param endness is the reads paired or unpaired
#' @param geneAnnotation reference annotation File
#' @return  txdb object to be used for alignment
#' @noRd
prepare_for_alignment<- function(outdir,filenames,SampleName,endness,entity,
                                geneAnnotation){
    trim.dir <- file.path(outdir , "fastp_results")
    writesampleFile(outdir, filenames,
                                SampleName,  endness)
    txdbfilename <- paste0(gsub(" ", "", entity), "_txdbobj", collapse = "")
    if (!file.exists(file.path(outdir, txdbfilename,
                            fsep = .Platform$file.sep))) {
        message("STEP 2; making txdb obj")
        txdb <-make_txdbobj(geneAnnotation, entity, outdir)
    } else {
        txdb <- AnnotationDbi::loadDb(file.path(outdir, txdbfilename,
                                                fsep = .Platform$file.sep))
    }
    return( txdb)
}

#' Run the alignment and count the genes
#' 
#' @param ref.dir : path to reference directory which contain
#' reference file(*.fa) and annotation file(*.gtf), can be NA
#' @param outdir : give the name of parent result dir , this can be existing or
#' not , rest of directory will be formed by program for organization
#' @param entity : Scientific name of species whose RNA is being analyzed
#' @param corenum : number of cores available for analysis #defaut 2
#' @param compare : what is sample/experimental design you have,
#'    paired or unpaired, as.group
#' @param aligner : One of "Rhisat2" or "Rbowtie2"; Rbowtie2 can be very slow
#'    for human and eukaryotic species
#' @param keep_tmp : set TRUE if keeping the aligned bam files, if set FALSE,
#'    bam files are deleted
#' @param cacheDir : directory where temporary files created during alignment
#' @param aligned_bam path where aligned bam is stored
#' @param endness is the reads PE or SE
#' @param sampleFile the path to sampleFile to use by qAlign
#' @param geneAnnotation geneAnnotation to use for counting genes
#' @param txdb txdb object for counting
#' @param genomeFile genomeFile for alignemnt
#' @return the analysis status message
#' @importFrom stringr str_replace_all
#' @param coldata phenofile imported in R
#' @import parallel
#' @import utils
#' @noRd
align_n_count <- function(aligned_bam,corenum, endness,keep_tmp, sampleFile, 
                                geneAnnotation, ref.dir, cacheDir, aligner,
                                outdir, txdb, entity, 
                                compare, genomeFile,coldata){
    if (!file.exists(file.path(aligned_bam, "alltrimmedalignedobj.RDS",
                            fsep = .Platform$file.sep))) {
        message("STEP 3 : aligning the sequence")
        aligned_proj <- run_qAlign(
            corenum, endness, sampleFile, genomeFile,
            geneAnnotation, ref.dir=ref.dir, cacheDir, aligner)
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
        ))) }
    
    if (compare == "paired"){
        coldata$SampleName <-  paste0("s", seq_len(length(coldata$SampleName))
                                    , "_", coldata$SampleName)
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
    return(list("counts" = as.data.frame(cnts),"grp.idx"= list(grp.idx)  ))
}

#' this function prepares and runs differential gene expression analysis if 
#' necessay
#' @param entity : Scientific name of species whose RNA is being analyzed
#' @param diff.tool : what differential tool to to use,
#' “DESEQ2” or “edgeR” available
#' @param compare : what is sample/experimental design you have,
#'    paired or unpaired, as.group
#' @param deseq2.dir directory path to store deseq2 results
#' @param cnts table of genes counts to do analysis
#' @param grp.idx the list of class of samples for comparision
#' @param SampleName name of sample to use for paired analysis
#' @param edger.dir where results of edgeR will be stored
#' @return exp.fc the fold change for each gene obtained from the analysis
#' @noRd
diff_analysis <- function(outdir,cnts, grp.idx,  entity, 
                        SampleName,compare,diff.tool ){
    deseq2.dir <- file.path(outdir, "differential_analysis/DESeq2")
    edger.dir <- file.path(outdir,"differential_analysis/edgeR")
    if (!file.exists(file.path(deseq2.dir, "Volcano_deseq2.tiff",
                            fsep = .Platform$file.sep))) {
        message("STEP 5a ; running differential analysis using DESeq2")
        exp.fcncnts.deseq2 <- run_deseq2(cnts, grp.idx, deseq2.dir, entity, 
                                        SampleName,compare)
    } else {
        deseq2.res.df <- read.table(
            file.path(deseq2.dir, "DESEQ2_logfoldchange.txt", 
                    fsep = .Platform$file.sep ),
            header = TRUE, sep = "\t", row.names = 1 )
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
    return( exp.fc)
}


#' This function runs different gene set analysis
#'
#' @param outdir : give the name of parent result dir , this can be existing or
#' not , rest of directory will be formed by program for organization
#' @param entity : Scientific name of species whose RNA is being analyzed

#' @param compare : what is sample/experimental design you have,
#'    paired or unpaired, as.group
#' @param cpd_id : compound id type
#' @param csamp  : sample column of compound
#' @param cref  : reference column of compound
#' @param mode :how to select pathways,"auto","combined"=using both cpd and gene
#' @param cdatapath path for the tsv file of compound
#' @param ccompare what is the comparision for sample and reference
#' @param pid pathway id
#' @param exp.fc logfold change of gene expression
#' @param cnts table of count of gene data with genes in row
#' @return returnval
#' @noRd
genecptsetanalysis <- function(outdir, entity,exp.fc, compare,  cnts,
                            cdatapath,cpd_id, csamp,cref,ccompare, mode,pid){
    gage.dir <- file.path(outdir, "gage_results")
    cset_dir <- file.path(gage.dir , "kegg_csets")
    combined_dir <- file.path(gage.dir, "combined_analysis_kegg")
    if (!file.exists(file.path(gage.dir,"gene_gageresults.rds", 
                            fsep = .Platform$file.sep))) {
        message("STEP 7 : running gene pathway analysis using GAGE")
        message(paste0(compare, "this is from pathwrap ", collapse = ""))
        res_gage_gene <-run_pathway(entity, exp.fc, compare, gage.dir, cnts)
        saveRDS(res_gage_gene, file.path(gage.dir,"gene_gageresults.rds"))
    } else {
        res_gage_gene <-readRDS(file.path(gage.dir,"gene_gageresults.rds"))
    }
    gpath_ids <- res_gage_gene$pathways_selected
    pgs_gene <- res_gage_gene$pgs.gene
    gage_out <- res_gage_gene$gage_result
    gsets <- res_gage_gene$gene_sets
    if (mode != "gene" ){
        if (mode == "compound_only" | mode == "combined"){
            if(!file.exists(file.path(cset_dir,"cpd_gageresults.rds",
                            fsep = .Platform$file.sep )) & !is.na(cdatapath)){
            message("STEP 8: running compound set analysis using GAGE")
            res_gage_cpd <-run_cpathway(cdatapath,cpd_id, csamp,cref,
                                    ccompare,cset_dir, entity )
        saveRDS(res_gage_cpd, file = file.path(cset_dir,"cpd_gageresults.rds"))
    } else {
        res_gage_cpd <-readRDS(file.path(cset_dir,"cpd_gageresults.rds"))
    }
    cpath_ids <- res_gage_cpd$pathways_selected
    pgs_cpd <- res_gage_cpd$pgs.gene
    gage_out_cpd <- res_gage_cpd$gage_result
    cpd_data <- res_gage_cpd$data_used
        }
    if (mode == "combined"){
        message("STEP 9: running combined gene set analysis using GAGE")
        qcut <- 0.01
        path_ids <- run_combinedpath_analysis(gpath_ids, cpath_ids,gsets,
                pgs_gene,pgs_cpd, combined_dir, gage_out, gage_out_cpd,qcut)
        plotpathways(combined_dir,entity,path_ids,
                    gene_data= exp.fc,cpd_data = cpd_data)
    }
    }
}    
