#' Runs fastp for quality and adapter trimming
#'
#' This function takes name of the samples and for each sample does the quality
#' and adapter trimming for Illumina and long read sequencing. It works for both
#' PE and SE data
#' @import ShortRead
#' @import Rfastp 
#' @import tools
#' @import reshape2
#' @import scales
#' @param outdir directory to store the trimmed fasta files
#' @param corenum number of threads 
#' @param phenofile : path to the phenofile where raw data files path is stored 
#' @return no value returned
#' @export
run_fastp <- function(phenofile, outdir, corenum) {
    #different adapters can be used for different seq_tech
    phenofile_res <- process_phenofile(phenofile)
    endness <- phenofile_res$endness
    FileName <- phenofile_res$FileName
    trim.dir <- file.path(outdir, "fastp_results")
    if (!file.exists(trim.dir)){ 
        file.create(trim.dir)
    }
    if (length(list.files(trim.dir, pattern = "_R1")) < length(phenofile_res$SampleName)){
        #adapterFasta <- ""
        df_res <- data.frame(row.names =c("passed_filter_reads", "low_quality_reads" ,  "too_many_N_reads" ,"too_short_reads",         
                                      "too_long_reads" , "adapter_trimmed_reads" , "adapter_trimmed_bases" ))
        for (indx in 1:length(phenofile_res$SampleName)){
            if (endness == "PE") 
            {
                trimmedoutfile1 <- file.path(trim.dir, 
                basename(file_path_sans_ext(file_path_sans_ext(
                file_path_sans_ext(FileName$FileName1[indx])))),
                                            fsep = .Platform$file.sep)
                trimmedoutfile2 <- file.path(trim.dir, 
                                             basename(file_path_sans_ext(file_path_sans_ext(
                                                 file_path_sans_ext(FileName$FileName2[indx])))),
                                             fsep = .Platform$file.sep)
                
                trimsumobj<-rfastp(path.expand(FileName$FileName1[indx]), 
                       path.expand(FileName$FileName2[indx]), path.expand(trimmedoutfile1),
                    thread = corenum, unpaired= file.path(trim.dir, "unpaired",
                            paste0(trimmedoutfile1, "unpaired.fastq.gz")) ,
                                   failedOut= file.path(trim.dir, "fastp_failed",
                                                        paste0(trimmedoutfile1, "fastp_failed.fastq.gz" )))
                
                file.rename(from = paste0(trimmedoutfile1, "_R2.fastq.gz"), 
                            to = paste0(trimmedoutfile2,"_R2.fastq.gz"))
                #need somthing here for plot of trim
                
            } else {
                 
                trimmedoutfile <- file.path(trim.dir, 
                    basename(file_path_sans_ext(file_path_sans_ext(FileName$FileName[indx]))),
                                            fsep = .Platform$file.sep)
                trimsumobj <- rfastp(path.expand(FileName$FileName[indx]), outputFastq = path.expand(trimmedoutfile),
                   thread =  corenum, failedOut= file.path(trim.dir, "fastp_failed",
                                                                          paste0(trimmedoutfile, "fastp_failed.fastq.gz" )))
                
            }
            trim_summary <- trimSummary(trimsumobj)[1:7,]
            df_res[phenofile_res$SampleName[indx]] <- as.numeric(trim_summary )
        }
        df_t <-as.data.frame(t(df_res))
        write.table(df_t, file = file.path(trim.dir, "trimming_result.txt") )
        plottrimresult(trim.dir, df_t)
    }

    return(invisible(NULL))
}
    
plottrimresult<- function(trim.dir,df_t ){
    df_t$sample <- rownames(df_t)
    df_long <- melt(df_t,id.vars="sample" , variable.name = "read_type", value.name = "count")
    
    # Compute proportion manually
    total_counts <- aggregate(count ~ sample, data = df_long, sum)
    df_long$proportion <- mapply(function(s, c) {
        total <- total_counts$count[total_counts$sample == s]
        c / total
    }, df_long$sample, df_long$count)
    
    
    tiff(file.path(trim.dir, "read_trimming_result_barplot.tiff",fsep = .Platform$file.sep),
         units = "in",  width = dim(df_t)[1], height = dim(df_t)[1],  res = 300)
    
    p<- ggplot(df_long, aes(x = sample, y = count, fill = read_type)) +
        geom_bar(stat = "identity", position = "fill") +
        geom_text(
            aes(label = ifelse(proportion > 0.05, percent(proportion, accuracy = 1), "")),
            position = position_fill(vjust = 0.5),
            size = 3,
            color = "white"
        ) +
        scale_y_continuous(labels = percent) +
        labs(y = "Proportion", x = "Sample", fill = "Read Type") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "bold"))+
        theme(axis.text.y = element_text(face = "bold"))
    plot(p)
    dev.off()
    return(invisible(NULL))
}
