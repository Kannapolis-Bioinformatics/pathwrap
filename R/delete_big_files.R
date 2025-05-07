
delete_tmp_files <- function(folder_of_tmp){
    ans <- readline(paste0("Are you sure you want to delete aligned reads ",
                           "yes " ,"or ",  "no "  ,  folder_of_tmp , "? "))
    if (substr(ans, 1, 1) == "n"){
        message("The aligned reads are not deleted")
    } else {
        unlink(list.files(folder_of_tmp, pattern = ".bam$|.fastq"
                          ,full.names = TRUE)) }
    
}