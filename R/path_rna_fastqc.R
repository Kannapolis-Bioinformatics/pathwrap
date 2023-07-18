#' Title
#'
#' @param fq.dir : the directory in which raw RNAseq files are stored
#' @param qc.dir : the directory in which results of quality check are stored
#' @param corenum : the number of cores used for running quality check
#'
#' @import ggplot2
#' @import fastqcr
#' @importFrom grDevices dev.off
#' @importFrom grDevices tiff
#' @export
#'

run_qc <- function(fq.dir, qc.dir, corenum){
  on.exit(closeAllConnections())
#install fastqc if system( "which fastqc", intern = T) fails
if (Sys.which("fastqc")=="" & !file.exists(paste0(qc.dir, "/FastQC/fastqc"))){
  ## work here
  fastqcr::fastqc_install(dest.dir = qc.dir)
  fastqc.path <- paste0(qc.dir, "/FastQC/fastqc")
}
else{
  fastqc.path <- Sys.which("fastqc")
}
  print("This is the fastqc tool we will run")
print(fastqc.path)
#print(paste0(fq.dir))
 #     , " ", qc.dir , " " ,unname(fastqc.path), "This is what i am checking " ))
fastqcr::fastqc(fq.dir, qc.dir,fastqc.path = unname(fastqc.path), threads = corenum) #what if threads is removed
print("Complete running fastqc")
qc <- qc_aggregate(qc.dir)
 print (" plotting total sequence and status of qc check in tiff files")
 #pdf(file.path(qc.dir,"total_seq.pdf"))# ,units="in", width=15, height=15, res=300)
 tiff( file.path(qc.dir,"total_seq.tiff") ,units="in", width=15, height=15, res=300)
 g <- ggplot(qc, aes(x=sample, y=tot.seq))+ geom_bar(stat="identity", position = "dodge")   + theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
 plot(g)
 dev.off()
# # pdf(file.path(qc.dir,"qc_heatmap.pdf"), width=15, height=15, res=300)
 tiff(file.path(qc.dir,"qc_heatmap.tiff"), units="in" , width=15, height=15, res=300)
 g <-    ggplot(qc, aes(x=module, y =sample, fill = status)) + geom_tile()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
 plot(g)
 dev.off()
 print("Multiqc report done")
}


