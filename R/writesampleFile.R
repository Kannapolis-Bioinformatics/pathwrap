#' Write sampleFile to run qAlign 
#'
#' This function creates the sampleFile to run qAlign
#' @param outdir :  path to result dir where the file is written
#' @param filenames : list of filenames extracted from phenofile
#' @return statement about the function run
writesampleFile <- function(outdir, filenames, SampleName, trim.dir, endness){
    sampleFile <- file.path(outdir, "sampleFile.txt")
    rawfileName <- as.data.frame(vapply(
        filenames, function(x) basename(x),
        character(dim(filenames)[1])
    ))
    fastp_files_name <- as.data.frame(vapply(
        as.list(rawfileName),
        function(x) str_replace_all(x, ".fastq.gz$", "_R1.fastq.gz"),
        character(dim(filenames)[1])
    ))
    FileName <- vapply(
        as.list(fastp_files_name),
        function(x) file.path(trim.dir, x),
        character(dim(filenames)[1])
    )
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
    return(sampleFile)
}
