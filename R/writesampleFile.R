#' Write sampleFile to run qAlign 
#'
#' This function creates the sampleFile to run qAlign
#' @import stringr
#' @param outdir :  path to result dir where the file is written
#' @param SampleName name of the sample to match to filenames
#' @param trim.dir directory to find trimmed reads
#' @param endness  paired or single end 
#' @param filenames : list of filenames extracted from phenofile
#' @return statement about the function run
writesampleFile <- function(outdir, filenames, SampleName, trim.dir, endness){
        sampleFile <- file.path(outdir, "sampleFile.txt")
        
        if (endness == "SE") {
            FileName <- file.path(trim.dir, 
                                    str_replace(SampleName,
                                            "$", "_R1.fastq.gz"))
            write.table(
                file = sampleFile, sep = "\t",
                as.data.frame(cbind(FileName, SampleName)),
                col.names = c("FileName", "SampleName"),
                quote = FALSE, row.names = FALSE
            )
        
    } else {
        FileName1 <- file.path(trim.dir, 
                            str_replace(SampleName,
                                        "$", "_R1.fastq.gz"))
        FileName2 <- file.path(trim.dir, 
                            str_replace(SampleName,
                            "$", "_R2.fastq.gz"))
        write.table(
            file = sampleFile, sep = "\t",
            as.data.frame(cbind(
                FileName1,
                FileName2, SampleName
            )),
            col.names = c("FileName1", "FileName2", "SampleName"),
            quote = FALSE, row.names = FALSE
        )
    }
    return(sampleFile)
}
