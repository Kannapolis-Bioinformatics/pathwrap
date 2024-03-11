# This code creates the phenofile and runs the pathviewrap

# install pathwrap from github or bioconductor and load library
# BiocManager::install("Kannapolis-Bioinformatics/pathwrap",
# force  = TRUE)
library(pathwrap)
# create directory to store results
Results <- tempdir()
phenofile <- tempfile()

# create columns for phenofile, this is for SE data
# Table col.names should be SampleName, FileName and Class
library(stringr)
FileName <- list.files(
    file.path(system.file(
    package = "pathwrap"
    ), "extdata"),
    full.names = TRUE,
    pattern = "fastq.gz"
)

patternmy <- c(dirname(FileName[1]), "_sub.fastq.gz")
SampleName <- str_replace_all(
    pattern = paste0(dirname(FileName)[1], "/|_sub.fastq.gz"),
    string = list.files(
    file.path(
      system.file(package = "pathwrap"), "extdata"
    ),
    full.names = TRUE, pattern = "fastq.gz"
    ), ""
)


Class <- c("A", "B", "A", "B")
write.table(as.data.frame(cbind(SampleName, FileName, Class)),
    file = phenofile, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE
)


message("this is the phenofile ", phenofile)
system.time({
    pathviewwrap(
    ref.dir = NA, phenofile = phenofile,
    outdir = Results, entity = "Mus musculus", corenum = 16,
    compare = "as.group",  keep_tmp = TRUE,
    startover = TRUE, diff.tool = "DESeq2", aligner = "Rhisat2"
    )
})
