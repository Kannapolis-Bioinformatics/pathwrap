#' Runs fastp for quality and adapter trimming
#'
#' This function takes name of the samples and for each sample does the quality
#' and adapter trimming for Illumina and long read sequencing. It works for both
#' PE and SE data
#' @import ShortRead
#' @import Rfastp 
#' @param FileName name of the raw files
#' @param endness wheatehr the raw data is paired or single ended
#' @param trim.dir directory to store the trimmed fasta files
#' @param corenum number of threads 
#' @param sampleName : name of the sample
#' @return no value returned

run_fastp <- function(sampleName, FileName, endness,
                    trim.dir, corenum) {
    #different adapters can be used for different seq_tech
    adapterFasta <- ""
    trimmedoutfile <- file.path(trim.dir, sampleName, fsep = .Platform$file.sep)
    infile <- FileName
    if (endness == "PE") {
        rfastp(path.expand(infile$FileName1), path.expand(infile$FileName2)
            , path.expand(trimmedoutfile),
            adapterFasta, thread = corenum)
    } else {
        rfastp(path.expand(infile), outputFastq = path.expand(trimmedoutfile),
            adapterFasta, thread =  corenum)
    }

    return(invisible(NULL))
}
