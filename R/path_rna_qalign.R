#############################################################################
#2. RUN THE ANALYSIS # Alignmnet and counting
#############################################################################

#' Title
#'
#' @param corenum : the number of cores used during alignment
#' @param endness : weather its paired end or single end
#' @param sampleFile : the file where information about location of sample is stored, see qAlign for more
#' @param genomeFile : the genome file used for alignment or bioconductor genome package
#' @param geneAnnotation : gene annnotaion file used for gene counting or bioconductor annotation package
#' @param ref.dir : directory in which genomeFile and genomeAnnotaion are stored
#' @param cacheDir : directory where temporary files genrated during alignment are store
#' @param  aligner : weather Rhisat2 or Rbowtie should be used for alignment
#'
#' @importFrom QuasR qAlign
#' @import Rhisat2
#'
#' @return aligned_proj : the R object generated from the alignment step
#' @export
#'

run_qAlign <- function(corenum, endness, sampleFile, genomeFile,geneAnnotation, ref.dir,cacheDir, aligner){
  #does ref.dir also have ref index, if not make indexes
  if( !is.na(ref.dir)){
      if( length(list.files(ref.dir , ".Rhisat2$", full.names = T)) !=1){
    sampleFiletmp <- read.table(sampleFile, "\t", header = T)[1,]
    sampleFiletmp_name <- paste0(gsub("sampleFile.txt","sampleFiletmp.txt", sampleFile))
    write.table(sampleFiletmp, sep=  "\t", col.names = T, row.names = F, file = sampleFiletmp_name)
    cl2 <- makeCluster(corenum)
    if(endness == "SE"){
      aligned_proj <-  QuasR::qAlign(sampleFiletmp_name, paired ="no", clObj=cl2, alignmentsDir =aligned_bam ,
                                     genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner =aligner,cacheDir= cacheDir)
    }else{
      aligned_proj <-  QuasR::qAlign(sampleFiletmp_name, paired ="fr", clObj=cl2, alignmentsDir =aligned_bam ,
                                     genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner =aligner ,cacheDir=cacheDir)
    }# this will form the reference index

    #the program check for aligned bam before running so we dont really need to remove this sample from our sampleFile
    unlink(sampleFiletmp_name)
    print("We made tmp file, and made index and one alignment")
      }
}
  cl2 <- makeCluster(corenum)
  print("Alignment is running")
  if (endness=="PE"){
    aligned_proj <- QuasR::qAlign(sampleFile,paired ="fr", clObj=cl2, alignmentsDir =aligned_bam , genome=genomeFile, geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner =aligner ,cacheDir=cacheDir)
  } else {

    aligned_proj <- QuasR::qAlign(sampleFile,paired ="no", clObj=cl2, alignmentsDir =aligned_bam ,  genome=genomeFile,geneAnnotation=geneAnnotation, splicedAlignment =TRUE, aligner =aligner,cacheDir =cacheDir)
  print("done")
    }
  saveRDS(aligned_proj, file.path(aligned_bam , "alltrimmedalignedobj.RDS"))
  stopCluster(cl2)
  return(aligned_proj)
}
