#' Run fastqc analysis in R
#'
#' runs fastqc analysis in R using fastqcr. If fastqc is not available in system
#' to run by R. this function is capable of downloading the fastqc tools before
#' running the quality check. The results of quality check is aggregated and
#' barplot of the total sequence and heatmap of the status of the qc check is
#' produced. do not install fastqc in your results directory
#'
#' @param fq.dir : the directory in which raw RNAseq files are stored
#' @param outdir : the directory in which results of quality check are stored
#' @param corenum : the number of cores used for running quality check
#' @return message of fastqc
#' @import dplyr
#' @import ggplot2
#' @import fastqcr
#' @importFrom grDevices dev.off
#' @importFrom grDevices tiff
#' @return message about the run if successfully
#' @export

run_qc <- function(fq.dir, outdir, corenum) {
    qc.dir <- file.path(outdir, "fastqc_results")
    if (!file.exists(file.path(qc.dir, "qc_heatmap.tiff",
                               fsep = .Platform$file.sep))){
        on.exit(closeAllConnections())
        # install fastqc if system( "which fastqc", intern = TRUE) fails
        
        if (Sys.which("fastqc") == "" & 
            !file.exists(file.path("~/bin", "FastQC", "fastqc",
                                fsep = .Platform$file.sep ))) {
            ## work here
            message("Please install fastqc using", "\n", 
                    "fastqcr::fastqc_install( )" , "\n", 
                    "Make sure it can be executed by R" )
            return(invisible(NULL))
        } else if (file.exists(file.path("~/bin", "FastQC", "fastqc",
                                        fsep = .Platform$file.sep))){
            fastqc.path <- "~/bin/FastQC/fastqc"
        } else {
            fastqc.path <- Sys.which("fastqc")
        }
    
        message("This is the fastqc tool we will run")
        message(fastqc.path)
        # check use of threadnum
        fastqcr::fastqc(fq.dir, qc.dir,
            fastqc.path = unname(fastqc.path), threads = corenum )
        message("Complete running fastqc")
        qc <- qc_aggregate(qc.dir)
        message(" plotting total sequence and status of qc check in tiff files")
    
        tiff(file.path(qc.dir, "total_seq.tiff", fsep = .Platform$file.sep),
            units = "in",width = dim(qc)[2], height = dim(qc)[2], res = 300 )
        
        g <- ggplot(qc, aes(x = sample, y = tot.seq)) +
            geom_bar(stat = "identity", position = "dodge", fill = "steelblue") +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,face = "bold"))+
            theme(axis.text.y = element_text(face = "bold"))
        plot(g)
        dev.off()
        # # pdf(file.path(qc.dir,"qc_heatmap.pdf"), width=15, height=15, res=300)
        tiff(file.path(qc.dir, "qc_heatmap.tiff",fsep = .Platform$file.sep),
            units = "in", width = 11,
            height = dim(qc)[2], res = 300
        )
        g <- ggplot(qc, aes(x = module, y = sample, fill = status)) +
            geom_tile() +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,size= 15,face = "bold"))+
            theme(axis.text.y = element_text(face = "bold"))
        plot(g)
        dev.off()
    
    
        all_fail_samples <- qc %>%
                 group_by(sample) %>%
                 summarize(all_fail = all(status == "FAIL")) %>%
                 filter(all_fail) %>%
             pull(sample)
        
        if (length(all_fail_samples)>0){
            return(invisible(NULL))
        }
    }
    return("Fastqc report done")
}
