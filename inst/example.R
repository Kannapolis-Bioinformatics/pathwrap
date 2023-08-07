# This code creates the phenofile and runs the pathviewrap

#install pathwrap from github or bioconductor and load library
#BiocManager::install("Kannapolis-Bioinformatics/pathwrap", force  = T)
library(pathviewwrap)
#create directory to store results
Results <- tempdir()
phenofile <-"hellotmpphenofile.txt"

#create columns for phenofile, this is for SE data
#Table col.names should be SampleName, FileName and Class
library(stringr)
FileName <- list.files(file.path(system.file(
    package = "pathviewwrap"), "extdata"), full.names = T)

patternmy <- c(dirname( FileName[1]) , "_sub.fastq.gz")
SampleName <- str_replace_all(
    pattern = paste0(dirname(FileName)[1], "/|_sub.fastq.gz" ),
    string =  list.files(file.path(
    system.file(package = "pathviewwrap"), "extdata"), 
    full.names = T) , "")


Class <- c("A", "B", "A", "B")
write.table(as.data.frame(cbind(SampleName, FileName, Class)), 
            file = phenofile, sep = "\t", row.names = F, 
            col.names = T, quote = F)


message("this is the phenofile ", phenofile )
system.time({
    pathviewwrap(
        ref.dir = NA, phenofile = phenofile,
        outdir = Results, entity = "Mus musculus", corenum = 16,
        compare = "as.group", seq_tech = "Illumina", keep_tmp = TRUE,
        rerun = FALSE, diff.tool = "DESeq2", aligner = "Rhisat2"
    )
})
