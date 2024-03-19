#' Write sampleFile to run qAlign 
#'
#' This function creates the sampleFile to run qAlign
#' @import stringr
#' @param outdir :  path to result dir where the file is written
#' @param SampleName name of the sample to match to filenames
#' @param trim.dir directory to find trimmed reads
#' @param endness  paired or single end 
#' @param filenames : list of filenames extracted from phenofile
#' @import tools
#' @return statement about the function run
writesampleFile <- function(outdir, filenames, SampleName, trim.dir, endness){
        sampleFile <- file.path(outdir, "sampleFile.txt", 
                                fsep = .Platform$file.sep)

        print("Above is str of filenames")
        if (endness == "SE") {
            print(filenames[,1])
            FileNametowrite <-  file.path(trim.dir , 
            str_replace_all(file_path_sans_ext(file_path_sans_ext(
                file_path_sans_ext(basename(filenames[,1])))), 
                        pattern = "$", replacement = "_R1.fastq.gz"))
            write.table(
                file = sampleFile, sep = "\t",
                as.data.frame(cbind(FileNametowrite, SampleName)),
                col.names = c("FileName", "SampleName"),
                quote = FALSE, row.names = FALSE
            )
        
    } else {
        FileName1 <- str_replace_all(file.path(trim.dir , 
            basename(file_path_sans_ext(file_path_sans_ext(file_path_sans_ext(
            filenames$FileName1)))), fsep = .Platform$file.sep), "$", 
            replacement = "_R1.fastq.gz")
        FileName2 <- str_replace_all(file.path(trim.dir , 
                    basename(file_path_sans_ext(file_path_sans_ext(file_path_sans_ext(
                    filenames$FileName2)))), fsep = .Platform$file.sep), "$", 
                    replacement = "_R2.fastq.gz")
        
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
