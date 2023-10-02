#' Runs fastp for quality and adapter trimming
#'
#' This function takes name of the samples and for each sample does the quality
#' and adapter trimming for Illumina and long read sequencing. It works for both
#' PE and SE data
#' @import ShortRead
#' @import Rfastp 
#' @param FileName name of the raw files
#' @param seq_tech long read or nanopore
#' @param endness wheatehr the raw data is paired or single ended
#' @param trim.dir directory to store the trimmed fasta files
#' @param corenum number of threads 
#' @param sampleName : name of the sample
#' @return no value returned

run_fastp <- function(sampleName, FileName, seq_tech, endness,
                    trim.dir, corenum) {
    # trimmedoutfile <- file.path(paste0(trim.dir, "/",
    #  sampleName, collapse = "|"))
    trimmedoutfile <- file.path(trim.dir, sampleName, fsep = .Platform$file.sep)
    infile <- FileName

    # actual command
    if (seq_tech == "PacBio" | seq_tech == "Nanopore") { # use custom adapters
        adapters_pacbio <- data(pacbioadapters,
            package = "pathviewwrap",
            envir = environment()
        )
        tmpadpfile <- tempfile("pacbioadapters.fna")
        writeFasta(file = tmpadpfile, adapters_pacbio)
        adapterFasta <- tmpadpfile
        adapterFasta <- system.file("extdata",
            package = "pathviewwrap",
            "adapters.fna"
        )
    } else {
        adapterFasta <- ""
    }
    message("STEP 1b : running fastp")
    if (endness == "PE") {
        #checkcondition <- !file.exists(trimmedoutfile["FileName1"]) &
            #!file.exists(trimmedoutfile["FileName2"])
        rfastp(path.expand(infile$FileName1), path.expand(infile$FileName2)
            , path.expand(trimmedoutfile),
            adapterFasta, thread = corenum)
    } else {
        #checkcondition <- !file.exists(trimmedoutfile)
        rfastp(path.expand(infile), outputFastq = path.expand(trimmedoutfile),
            adapterFasta, thread =  corenum)
    }

    return(invisible(NULL))
}
