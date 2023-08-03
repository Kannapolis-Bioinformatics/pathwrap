#' Runs fastp for quality and adapter trimming
#'
#' This function takes name of the samples and for each sample does the quality
#' and adapter trimming for Illumina and long read sequencing. It works for both
#'  PE and SE data
#'
#' @param sampleName : name of the sample
#' @return no value returned

run_fastp <- function(sampleName) {
    # integer
    intformatch <- grep(sampleName, FileName[, 1], value = FALSE, fixed = TRUE)
    # trimmedoutfile$FileName1 , trimmedoutfile$Filename2
    # infile$FileName1 , infile$FileName2 #infile for SE
    trimmedoutfile <- FileName[intformatch, ]
    infile <- filenames[intformatch, ]
    if (endness == "PE") {
        infileoutfile <- paste0(
            "-i ", infile$FileName1, " -I ",
            infile$FileName2, " -o ",
            trimmedoutfile["FileName1"], " -O ",
            trimmedoutfile["FileName2"]
        )
    } else {
        infileoutfile <- paste0("-i ", infile, " -o ", trimmedoutfile)
    }
    logfile <- paste0(
        " -h ", file.path(trim.dir, paste0(sampleName, ".html")),
        " -j ", file.path(trim.dir, paste0(sampleName, ".json"))
    )

    # actual command
    if (seq_tech == "PacBio" | seq_tech == "Nanopore") { # use custom adapters
        cmd <- paste0(
            "fastp ", infileoutfile, "--adapter_fasta",
            "data/adapters.fna", logfile
        )
    } else {
        cmd <- paste0("fastp ", infileoutfile, logfile)
    }
    # file check before running command
    # if (endness=="PE"){
    # checkcondition <- length(list.files(trim.dir, pattern ="json")) <=0 |
    # (length(list.files(trim.dir, pattern ="trimmed")) !=
    # (2*length(list.files(trim.dir, pattern ="json"))))
    # } else { checkcondition <- length(list.files(trim.dir, pattern ="json"))
    #<=0 | (length(list.files(trim.dir, pattern ="trimmed")) !=
    # length(list.files(trim.dir, pattern ="json")))}
    # if (checkcondition){
    #  message("STEP 1b : running fastp")
    # message(cmd)
    # system(cmd)
    if (endness == "PE") {
        checkcondition <- !file.exists(trimmedoutfile["FileName1"]) &
            !file.exists(trimmedoutfile["FileName2"])
    } else {
        checkcondition <- !file.exists(trimmedoutfile)
    }
    if (checkcondition) {
        message("STEP 1b : running fastp")
        message(cmd)
        system(cmd)
    }
    return(invisible(NULL))
}
