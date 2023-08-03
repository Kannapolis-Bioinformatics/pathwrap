#' 2.  RUN THE ANALYSIS # Alignment and counting
#'
#' this function runs the Rhisat2 or Rbowtie for alignment on paired or
#' single end mode. It saves the alignment object in RDS file which can
#' be loaded in R for further analysis. If the reference index is not found
#' in the reference directory, it creates reference index before running
#' alignment. If the reference genome is a package the reference index is
#' created as R package.It generates the barplot of mapped and unmapped sequence
#' reads.
#'
#' @param corenum : the number of cores used during alignment
#' @param endness : weather its paired end or single end
#' @param sampleFile : the file where information about location of sample is
#' stored, see qAlign for more
#' @param genomeFile : the genome file used for alignment or bioconductor genome
#'  package
#' @param geneAnnotation : gene annnotaion file used for gene counting or
#' bioconductor annotation package
#' @param ref.dir : directory in which genomeFile and genomeAnnotaion are stored
#' @param cacheDir : directory where temporary files generated during alignment
#' are store
#' @param  aligner : weather Rhisat2 or Rbowtie should be used for alignment
#'
#' @import QuasR
#' @importFrom QuasR qAlign
#' @import Rhisat2
#' @importFrom grDevices tiff
#' @import ggplot2
#'
#' @return R object generated from the alignment step
#'

run_qAlign <- function(corenum, endness, sampleFile, genomeFile, geneAnnotation,
                    ref.dir, cacheDir, aligner) {
    # does ref.dir also have ref index, if not make indexes
    if (!is.na(ref.dir)) {
        if (length(list.files(ref.dir, ".Rhisat2$", full.names = TRUE)) != 1) {
            sampleFiletmp <- read.table(sampleFile, "\t", header = TRUE)[1, ]
            sampleFiletmp_name <- paste0(gsub(
                "sampleFile.txt",
                "sampleFiletmp.txt", sampleFile
            ))
            write.table(sampleFiletmp,
                sep = "\t", col.names = TRUE,
                row.names = FALSE, file = sampleFiletmp_name
            )
            cl2 <- makeCluster(corenum)
            if (endness == "SE") {
                aligned_proj <- qAlign(sampleFiletmp_name,
                    paired = "no", clObj = cl2, alignmentsDir = aligned_bam,
                    genome = genomeFile, geneAnnotation = geneAnnotation,
                    splicedAlignment = TRUE, aligner = aligner,
                    cacheDir = cacheDir
                )
            } else {
                aligned_proj <- qAlign(sampleFiletmp_name,
                    paired = "fr", clObj = cl2, alignmentsDir = aligned_bam,
                    genome = genomeFile, geneAnnotation = geneAnnotation,
                    splicedAlignment = TRUE, aligner = aligner,
                    cacheDir = cacheDir
                )
            } # this will form the reference index

            # the program check for aligned bam before running so we dont really
            # need to remove this sample from our sampleFile
            unlink(sampleFiletmp_name)
            message("We made tmp file, and made index and one alignment")
        }
    }
    cl2 <- makeCluster(corenum)
    message("Alignment is running")
    if (endness == "PE") {
        aligned_proj <- qAlign(sampleFile,
            paired = "fr", clObj = cl2,
            alignmentsDir = aligned_bam, genome = genomeFile,
            geneAnnotation = geneAnnotation,
            splicedAlignment = TRUE,
            aligner = aligner, cacheDir = cacheDir
        )
    } else {
        aligned_proj <- qAlign(sampleFile,
            paired = "no", clObj = cl2,
            alignmentsDir = aligned_bam, genome = genomeFile,
            geneAnnotation = geneAnnotation,
            splicedAlignment = TRUE,
            aligner = aligner, cacheDir = cacheDir
        )
        message("done")
    }
    saveRDS(aligned_proj, file.path(aligned_bam, "alltrimmedalignedobj.RDS"))
    stopCluster(cl2)
    # Plot the alignment mapping statistics
    aligned_stat_my <- alignmentStats(aligned_proj)
    typesofdata <- c(
        rep("mapped", dim(aligned_stat_my)[1]),
        rep("unmapped", dim(aligned_stat_my)[1])
    )
    genomeofsamples <- c(rep(rownames(aligned_stat_my), 2))
    value <- c(aligned_stat_my[, 2], aligned_stat_my[, 3])
    data <- data.frame(genomeofsamples, typesofdata, value)
    tiff(file.path(aligned_bam, "mapping_stats.tiff"),
        units = "in", width = 15, height = 15, res = 300
    )
    g <- ggplot(data, aes(
        fill = typesofdata, y = value,
        x = stringr::str_remove_all(genomeofsamples, ":genome")
    )) +
        geom_bar(position = "fill", stat = "identity") +
        ylab("Proportion") +
        xlab("samples") +
        theme(legend.title = element_blank()) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    plot(g)
    dev.off()
    return(aligned_proj)
}
