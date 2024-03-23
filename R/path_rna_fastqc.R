#' Run fastqc analysis in R
#'
#' runs fastqc analysis in R using fastqcr. If fastqc is not available in system
#' to run by R. this function is capable of downloading the fastqc tools before
#' running the quality check. The results of quality check is aggregated and
#' barplot of the total sequence and heatmap of the status of the qc check is
#' produced. do not install fastqc in your results directory
#'
#' @param fq.dir : the directory in which raw RNAseq files are stored
#' @param qc.dir : the directory in which results of quality check are stored
#' @param corenum : the number of cores used for running quality check
#' @return message of fastqc
#'
#' @import ggplot2
#' @import fastqcr
#' @importFrom grDevices dev.off
#' @importFrom grDevices tiff
#'

run_qc <- function(fq.dir, qc.dir, corenum) {
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
        units = "in",
        width = 15, height = 15, res = 300
    )
    g <- ggplot(qc, aes(x = sample, y = tot.seq)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    plot(g)
    dev.off()
    # # pdf(file.path(qc.dir,"qc_heatmap.pdf"), width=15, height=15, res=300)
    tiff(file.path(qc.dir, "qc_heatmap.tiff",fsep = .Platform$file.sep),
        units = "in", width = 15,
        height = 15, res = 300
    )
    g <- ggplot(qc, aes(x = module, y = sample, fill = status)) +
        geom_tile() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    plot(g)
    dev.off()
    return("Multiqc report done")
}
