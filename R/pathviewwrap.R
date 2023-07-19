#' Title
#'
#' @param ref.dir : path to reference directory which contain reference file(*.fa) and annotation file(*.gtf), can be NULL
#' @param phenofile : path to phenofile ; see note on test data to see the format of phenofile
#' @param outdir  : give the name of parent result dir , this can be existing or not , rest of directory will be formed by program for organization
#' @param entity  : Scientific name of species whose RNA is being analyzed
#' @param corenum : number of cores available for analysis #defaut 2
#' @param diff.tool :  what differential tool to to use, “DESEQ2” or “edgeR” available
#' @param compare : what is sample/experimental design you have, paired or unpaired, as.group
#' @param seq_tech : Illumina, pacbio or nanopore
#' @param aligner : One of "Rhisat2" or "Rbowtie2"; Rbowtie2 can be very slow for human and eukaryotic species
#' @param keep_tmp : set TRUE if keeping the aligned bam files, if set FALSE, bam files are deleted
#' @param rerun : if FALSE the previously complete step will not be rerunned, if TRUE analysis starts from first step
#' @param cacheDir : directory where temporary files created during alignment are stored
#'
#'
#' @importFrom stringr str_replace_all
#' @import parallel
#' @import utils
#'
#' @export
pathviewwrap <- function( ref.dir = NA, phenofile= NA, outdir="results",  entity="Mus musculus",
                         corenum = 2,  compare="unpaired",diff.tool = "DESeq2", seq_tech="Illumina", keep_tmp = FALSE,rerun = FALSE, cacheDir = NULL , aligner){

  on.exit(closeAllConnections())
  aligned_bam <- NA
   dirlist <- sanity_check( ref.dir , pos= 1, outdir,  entity , corenum , compare, rerun)
   qc.dir <- dirlist[1]
   trim.dir <- dirlist[2]
   genomeFile <- dirlist[3]
   print("this is genome File")
   print(genomeFile)
   geneAnnotation <- dirlist[4]
   print("this is geneAnnotation")
   print(geneAnnotation)
  deseq2.dir <- dirlist[5]
  edger.dir <- dirlist[6]
  gage.dir <- dirlist[7]

  #duplicated codes
  # if (!file.exists(phenofile)){ ###TO DO make sure reference is first aplhanumerically#
  #   print("Please provide phenofile with Class information")
  # }
  # coldata <- read.table(phenofile, sep = "\t", header = T)
  # if(colnames(coldata)[ncol(coldata)]!="Class"){
  #   print("Please make sure class information is in last column with colname 'Class' . ")
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

  #run the fastqc


  if (!file.exists(phenofile)){ ###TO DO make sure reference is first aplhanumerically#
    print("Please provide phenofile with Class information")
  }
  coldata <- read.table(phenofile, sep = "\t", header = T)
  if(colnames(coldata)[ncol(coldata)]!="Class"){
    print("Please make sure class information is in last column with colname 'Class' . ")
  }
  coldata$Class <- as.factor(coldata$Class)
  SampleName <- coldata$Sample
  filenames <- as.data.frame(coldata[,-c(1,ncol(coldata))])

  if(dim(filenames)[2] == 1){
    endness <- "SE"
    fq.dir <-  dirname(filenames[1,1])
  } else if(dim(filenames)[2] == 2){
    endness <- "PE"
    fq.dir <- dirname(filenames$FileName1[1])
  }

  if (!file.exists(file.path(qc.dir,"qc_heatmap.tiff"))){
    print("STEP 1 ; running fastqc")
    print("this is qc.dir")
    print(qc.dir)
    run_qc(fq.dir, qc.dir, corenum)
  }

  sampleFile <- file.path(outdir, "sampleFile.txt")
  rawfileName <- as.data.frame(sapply(filenames, function(x) basename(x)))
  fastp_files_name <- as.data.frame(sapply(rawfileName, function(x) str_replace_all(x, ".fastq.gz$" ,"_trimmed.fastq.gz" )))
  FileName <- sapply(fastp_files_name, function(x) file.path(trim.dir , x))

  if(endness == "SE") {
    write.table(file =sampleFile,sep = "\t", as.data.frame( cbind( FileName, SampleName)) ,  col.names = c("FileName", "SampleName"),quote =F ,row.names=F)
  } else{
    write.table(file=sampleFile,sep = "\t", as.data.frame( cbind(FileName[,1], FileName[,2], SampleName)), col.names = c("FileName1","FileName2", "SampleName"), quote =F ,  row.names=F)
  }

  #just in case there is random component in run_fastp
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)

  cl <- makeCluster(corenum)
  clusterExport(cl,c("seq_tech","endness", "FileName" ,"filenames", "trim.dir"), envir = environment())#.GlobalEnv) ??
  parSapply(cl , SampleName  ,run_fastp )
  print("the trim run is complete")
  stopCluster(cl)

  #to check if all the nodes run fine
 # bad <- sapply(r, inherits, what = "try-error") # r<- mclappy()

  #make txdb from annotation
  if(!file.exists(paste0(outdir,"/", gsub(" ", "", entity), "_txdbobj"))){
    print("STEP 2; making txdb obj")
    txdb <- make_txdbobj(geneAnnotation, corenum, genomeFile, entity, outdir)
  } else{
    txdb <- AnnotationDbi::loadDb(paste0(outdir,"/", gsub(" ", "", entity), "_txdbobj"))
  }

  if(!file.exists(file.path(aligned_bam , "alltrimmedalignedobj.RDS"))){
    setwd(outdir)#make sure you delete this file before rerunning can be better
    print("STEP 3 : aligning the sequence")
    aligned_proj <- run_qAlign(corenum, endness, sampleFile, genomeFile,geneAnnotation, ref.dir, cacheDir,aligner) #can be better??
  }
  else{
    aligned_proj <- readRDS(file.path(aligned_bam , "alltrimmedalignedobj.RDS"))
  }

  if(!file.exists(file.path(outdir, "combinedcount.trimmed.RDS")  ))  {
    print("STEP 4: counting aligned sequences")
    cnts <-run_qCount( aligned_proj, corenum, outdir, txdb, entity)
    } else{
    cnts <- as.data.frame(readRDS(file.path(outdir, "combinedcount.trimmed.RDS") ))
    }
  #print("these are samplename for cnts , cnts[, coldata$SampleName] ")
  #print(coldata$SampleName)
  cnts <- cnts[, coldata$SampleName]
  if(  all(coldata$SampleName == colnames(cnts)) ){#if this then proceed
    ref <- which(coldata$Class ==  levels(as.factor(coldata$Class))[1])
    samp <- which(coldata$Class ==  levels(as.factor(coldata$Class))[2])
    grp.idx <-NULL
    grp.idx[ref] <- "reference"
    grp.idx[samp] <- "sample"

  } else{
    message("make sure pheno file have only samples analysed")
  }


  if (keep_tmp == FALSE){
    print("deleting aligned bam files, bam file index and log files")
    #unlink(file.path(outdir, "aligned_bam", "*bam*"))
    unlink(list.files(file.path(outdir, "aligned_bam"), pattern = ".bam$|.bai$", full.names = T))
  }

  if(!file.exists(paste0(deseq2.dir, "/Volcano_deseq2.tiff"))){
    print("STEP 5a ; running differential analysis using DESeq2")
    exp.fcncnts.deseq2 <- run_deseq2(cnts,grp.idx, deseq2.dir)
  }    else{
    deseq2.res.df  <- read.table(file.path(deseq2.dir, "DESEQ2_logfoldchange.txt"), header = T, sep = "\t", row.names = 1) #works with gage
    exp.fcncnts.deseq2 <- deseq2.res.df  $log2FoldChange
    names( exp.fcncnts.deseq2) <-  rownames(deseq2.res.df )

  }
  if(!file.exists(paste0(edger.dir, "Volcano_edgeR.tiff"))){
    print("STEP 5b ; running differential analysis using edgeR")
    exp.fcncnts.edger <- run_edgeR(cnts,grp.idx, edger.dir)
  } else{
    edger.res.df  <- read.table(file.path(edger.dir, "edgeR_logfoldchange.txt"), header = T, sep = "\t", row.names = 1) #works with gage
    exp.fcncnts.deseq2 <-edger.res.df  $log2FC
    names( exp.fcncnts.deseq2) <-  rownames(edger.res.df )
  }

  setwd(gage.dir)
  #chosing to use deseq2 result or edger result for gage
  if(diff.tool == "DESeq2"){
    exp.fc <- exp.fcncnts.deseq2
  } else{
    exp.fc <- exp.fcncnts.edger
  }
  if(!file.exists("*.txt")){
    print("STEP 6 : running pathway analysis using GAGE")
    print(paste0(compare, "this is from pathviewwrap"))
    run_pathway(entity,exp.fc, compare, gage.dir, cnts, grp.idx)
  }
}

#when loading pathview??
#####Loading required namespace: org.Mm.eg.db
#Installing package(s) 'org.Mm.eg.db'
#installing the source package ‘org.Mm.eg.db’
