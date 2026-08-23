#' Write sampleFile to run qAlign 
#'
#' This function creates the sampleFile to run qAlign
#' @import stringr
#' @param outdir path to result dir where the file is written
#' @param phenofile_res : object of process phenofile
#' @param nchunks how many pieces should the phenofile be cut to
#' @import tools
#' @import stringr
#' @return sampleFile sampleFile path to be used by qAlign

writesampleFile <- function(outdir, phenofile_res,nchunks){
    trim.dir <- file.path(outdir , "fastp_results")
    
    #phenofile_res <- process_phenofile(phenofile)
    if(nchunks>1){
        chunk_indices <- cut(seq_along(phenofile_res$SampleName), breaks = nchunks, labels = FALSE)
        index_chunks <- split(seq_along(phenofile_res$SampleName), chunk_indices)
        message(index_chunks)
    } else{
        index_chunks<- list(seq_along(phenofile_res$SampleName))
    }
    sampleFile_lst<- c()
    for (val in 1:length(index_chunks)){
        file.remove(list.files(tempdir()), recursive = TRUE)
        idxchunk<- index_chunks[[val]]
        message(idxchunk)
        sampleFile <- file.path(outdir, paste0( "sampleFile_", idxchunk[1] , "_", idxchunk[length(idxchunk)], ".txt"), 
                                fsep = .Platform$file.sep)
        sampleFile_lst<- append(sampleFile_lst,sampleFile)
        
        if (is.null(phenofile_res$coldata$Class)){
            sampleFile <- phenofile
        } else {
            SampleName <- phenofile_res$SampleName[idxchunk]
            filenames <- phenofile_res$FileName[idxchunk, , drop = FALSE]
            if (phenofile_res$endness == "SE") {
                FileNametowrite <-  file.path(trim.dir ,
                                              str_replace_all(file_path_sans_ext(file_path_sans_ext(
                                                  file_path_sans_ext(basename(filenames$FileName)))),
                                                  pattern = "$", replacement = "_R1.fastq.gz"))
              
                write.table(file = sampleFile, sep = "\t",
                            as.data.frame(cbind(as.data.frame(FileNametowrite), 
                                                SampleName)),col.names = c("FileName", "SampleName"),
                            quote = FALSE, row.names = FALSE )
            } else {
                
                FileName1 <- str_replace_all(file.path(trim.dir , 
                                                       basename(file_path_sans_ext(file_path_sans_ext(file_path_sans_ext(
                                                           filenames$FileName1)))), fsep = .Platform$file.sep), "$", 
                                             replacement = "_R1.fastq.gz")
                FileName2 <- str_replace_all(file.path(trim.dir , 
                                                       basename(file_path_sans_ext(file_path_sans_ext(
                                                           file_path_sans_ext(filenames$FileName2)))), 
                                                       fsep = .Platform$file.sep),"$",replacement = "_R2.fastq.gz")
                
                write.table(
                    file = sampleFile, sep = "\t",
                    as.data.frame(cbind(
                        FileName1,
                        FileName2, SampleName
                    )),
                    col.names =  c("FileName1","FileName2", "SampleName")
                    , quote = FALSE, row.names = FALSE        )
            }
        }
        
    }
    
    return(sampleFile_lst)
}

