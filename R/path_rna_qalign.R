#' 2. RUN THE ANALYSIS # Alignment and counting
#'
#' this function runs the Rhisat2 or Rbowtie for alignment on paired or
#' single end mode. It saves the alignment object in RDS file which can
#' be loaded in R for further analysis. If the reference index is not found
#' in the reference directory, it creates reference index before running
#' alignment. If the reference genome is a package the reference index is
#' created as R package.It generates the barplot for number of mapped and 
#' unmapped sequence reads.
#' @param phenofile_res : obje of process phenofile
#' @param corenum : the number of cores used during alignment
#' @param nchunks the number of segment you want to run the alignemnt
#' @param references : the list of references genome file and gene ananotaion used for 
#' alignment or bioconductor genome package
#' @param cacheDir : directory where temporary files generated during alignment
#' are store
#' @param aligner : weather Rhisat2 or Rbowtie should be used for alignment
#' @param outdir directory where aligned bam files are stored
#' @import QuasR
#' @importFrom QuasR qAlign
#' @import scales
#' @import Rhisat2
#' @importFrom grDevices tiff
#' @import ggplot2
#' @import dplyr
#' @return R object generated from the alignment step
#'

run_qAlign <-function(phenofile_res , cacheDir, aligner, references, outdir,corenum,nchunks){
    # does ref.dir also have ref index, if not make indexes
    aligned_bam<- file.path(outdir, "aligned_bam")
    setwd(aligned_bam)
    sampleFile_lst <- writesampleFile(outdir, phenofile_res,nchunks)
    
    genomeFile <- references$genomeFile
    geneAnnotation <- references$geneAnnotation
    
    if (dim(read.table(sampleFile_lst[1], "\t", header = T))[2]==3){
        pairedinfo <- "fr"
    } else{
        pairedinfo <- "no"
    }
    if (!file.exists(file.path(aligned_bam, "alltrimmedalignedobj.RDS",
                               fsep = .Platform$file.sep))){
        unlink(list.files(cacheDir, full.names = TRUE), recursive = T)
        aligned_proj_list <- c()
        data_all<- c()
        for (sampleFile in sampleFile_lst){
    
            cl2 <- makeCluster(corenum)
            message("Alignment is running")
    
            aligned_proj <- qAlign(sampleFile,
                paired = pairedinfo, clObj = cl2,
                alignmentsDir = aligned_bam, genome = genomeFile,
                geneAnnotation = geneAnnotation,
                splicedAlignment = TRUE,
                aligner = aligner, cacheDir = cacheDir)
            stopCluster(cl2)
            aligned_proj_list<-BiocGenerics::append( aligned_proj_list,aligned_proj)
            aligned_stat_my <- alignmentStats(aligned_proj)
            typesofdata <- c(
                rep("mapped", dim(aligned_stat_my)[1]),
                rep("unmapped", dim(aligned_stat_my)[1])
            )
            genomeofsamples <- c(rep(rownames(aligned_stat_my), 2))
            value <- c(aligned_stat_my[, 2], aligned_stat_my[, 3])
            data <- data.frame(genomeofsamples, typesofdata, value)
        
            data_all <- rbind(data_all, data)
            
        }
        print(data_all)
        plotalignmentstats(data_all, aligned_bam)
        saveRDS(aligned_proj_list, file.path(aligned_bam, "alltrimmedalignedobj.RDS",
                                             fsep = .Platform$file.sep))
        
    } else {
        aligned_proj_list<-  readRDS(  file.path(aligned_bam, 
                                             "alltrimmedalignedobj.RDS",
                                             fsep = .Platform$file.sep))
    }
    return(aligned_proj_list)
}

#' Plot the alignment mapping statistics
#' @param data_to_plot data to plot alignment stats
#' @param aligned_bam directory for bam files
#' @return  returns after completing plot
plotalignmentstats <- function(data_to_plot, aligned_bam){
    if (!file.exists(file.path(aligned_bam, "mapping_stats.tiff",fsep = .Platform$file.sep))){
        data_prop <- data_to_plot %>%
            group_by(genomeofsamples) %>%
            mutate(prop = value / sum(value))
        
        tiff(file.path(aligned_bam, "mapping_stats.tiff",fsep = .Platform$file.sep),
            units = "in", width = length(unique(data_to_plot$genomeofsamples)), height = 5, res = 300
        )
        g <- ggplot(data_to_plot, aes(
            fill = typesofdata, y = value,
            x = stringr::str_remove_all(genomeofsamples, ":genome")
        )) +
            geom_bar(position = "fill", stat = "identity") +
            geom_text(data=data_prop,
                aes(label = percent(prop, accuracy = 1)),
                position = position_fill(vjust = 0.5),
                size = 3,
                color = "steelblue"
            ) +
            ylab("Proportion") +
            xlab("samples") +
            theme(legend.title = element_blank()) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        print(g)
        dev.off()
    }
        return(invisible(NULL)) 
}

#use this if making for hpc cluster
# ref.dir <- dirname(references$genomeFile)
# if (!is.na(ref.dir)) {
#     apattern <- paste0(".", aligner, "$")
#     if (length(list.files(ref.dir, apattern, full.names = TRUE)) != 1) {
#         sampleFiletmp <- read.table(sampleFile, "\t", header = TRUE)[1, ]
#         sampleFiletmp_name <- paste0(gsub(
#             "sampleFile.txt",
#             "sampleFiletmp.txt", sampleFile))
#         write.table(sampleFiletmp,
#             sep = "\t", col.names = TRUE,
#             row.names = FALSE, file = sampleFiletmp_name)
#         cl2 <- makeCluster(corenum)
#         if (endness == "PE") {
#             pairedinfo <- "fr"
#         } else{ 
#             pairedinfo <- "no" }
#         aligned_proj <- qAlign(sampleFiletmp_name,
#         paired = pairedinfo, clObj = cl2,alignmentsDir = aligned_bam, 
#         genome = genomeFile, geneAnnotation = geneAnnotation,
#         splicedAlignment = TRUE,aligner = aligner, cacheDir = cacheDir) 
#         # this will form the reference index
#         # the program check for aligned bam before running so we dont really
#         # need to remove this sample from our sampleFile
#         unlink(sampleFiletmp_name)
#         message("We made tmp file, and made index and one alignment")
#     }
# }
