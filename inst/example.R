# This code creates the phenofile and runs the pathwrap
#Make sure required packages are installed if reference directory is not specified
#library(pathwrap)
#data(anntpkglist)
#genomepkg <- anntpkglist$genome[which(anntpkglist$species=="Mus musculus")]
#anntpkg <- anntpkglist$annotation[which(anntpkglist$species=="Mus musculus")]
#BiocManager::install(genomepkg)
#BiocManager::install(anntpkg)
# install pathwrap from github or bioconductor and load library
library(pathwrap)
# create directory to store results
#Results <- "/Users/edhungel/Documents/Research/myresults"
#phenofile <- "/Users/edhungel/Research/Documents/myphenofile.txt"
#provide actual path like above
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
csamp <- c(1,2)
cref <- c(3,4)
cpath <- file.path(system.file(package = "pathwrap", "extdata"), 
                   "example_cpd_data.tsv")
                   
message("this is the phenofile ", phenofile)
if(interactive()){ system.time({
  pathwrap(
    ref.dir = NA, phenofile = phenofile,mode = "gene", 
    outdir = Results, entity = "Mus musculus", corenum = 16,
    compare = "as.group",  keep_tmp = TRUE,
    startover = TRUE, diff.tool = "DESeq2", aligner = "Rhisat2",
    cdatapath = cpath, cref= cref, csamp = csamp, ccompare = "paired"
  )
})}
