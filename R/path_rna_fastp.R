#' Runs fastp for quality and adapter trimming
#'
#' This function takes name of the samples and for each sample does the quality
#' and adapter trimming for Illumina and long read sequencing. It works for both
#' PE and SE data
#' @import ShortRead
#' @import Rfastp 
#' @import tools
#' @param FileName name of the raw files, dataframe with colnames 
#' FileName1 and FileName2 for paired read
#' @param endness wheatehr the raw data is paired or single ended
#' @param trim.dir directory to store the trimmed fasta files
#' @param corenum number of threads 
#' @param sampleName : name of the sampleverbiage
#' @return no value returned

run_fastp <- function(sampleName, FileName, endness,
                    trim.dir, corenum) {
    #different adapters can be used for different seq_tech
    adapterFasta <- ""
    if (endness == "PE") {
        trimmedoutfile <- file.path(trim.dir, 
        basename(file_path_sans_ext(file_path_sans_ext(
        file_path_sans_ext(FileName$FileName1)))),
                                    fsep = .Platform$file.sep)
        rfastp(path.expand(FileName$FileName1), path.expand(FileName$FileName2)
            , path.expand(trimmedoutfile),
            adapterFasta, thread = corenum)
    } else {
        trimmedoutfile <- file.path(trim.dir, 
            basename(file_path_sans_ext(file_path_sans_ext(FileName))),
                                    fsep = .Platform$file.sep)
        rfastp(path.expand(FileName), outputFastq = path.expand(trimmedoutfile),
            adapterFasta, thread =  corenum)
    }

    return(invisible(NULL))
}
