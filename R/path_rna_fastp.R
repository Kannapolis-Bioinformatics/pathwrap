#' Runs fastp for quality and adapter trimming
#'
#' This function takes name of the samples and for each sample does the quality
#' and adapter trimming for Illumina and long read sequencing. It works for both
#' PE and SE data
#' @import ShortRead
#' @import Rfastp
#' @param sampleName : name of the sample
#' @return no value returned

run_fastp <- function(sampleName, FileName, seq_tech, endness) {
    print(str(FileName))
    # integer
    intformatch <- grep(sampleName, FileName[, 1], value = FALSE, fixed = TRUE)
    # trimmedoutfile$FileName1 , trimmedoutfile$Filename2
    # infile$FileName1 , infile$FileName2 #infile for SE
    trimmedoutfile <- str_remove(FileName[intformatch, 1], ".fastq.gz")
    infile <- FileName[intformatch, ]


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
        checkcondition <- !file.exists(trimmedoutfile["FileName1"]) &
            !file.exists(trimmedoutfile["FileName2"])
        rfastp(infile$FileName1, infile$FileName2, trimmedoutfile, adapterFasta)
    } else {
        checkcondition <- !file.exists(trimmedoutfile)
        rfastp(infile, outputFastq = trimmedoutfile, adapterFasta)
    }

    return(invisible(NULL))
}
