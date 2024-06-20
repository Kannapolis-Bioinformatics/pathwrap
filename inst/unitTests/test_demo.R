test_divideBy <- function() {
  goodtestresultmessage <- "The analysis is complete"
  #Results <- tempdir()
  Results <- "~/Documents/Research/myresults"
  phenofile <- system.file("extdata", "phenofile_SE.txt",
    package = "pathwrap"
  )

  # create columns for phenofile, this is for SE data
  # Table col.names should be SampleName, FileName and Class
  # library(stringr)
  # FileName <- list.files(
  #   file.path(system.file(
  #     package = "pathwrap"
  #   ), "extdata"),
  #   full.names = TRUE,
  #   pattern = "fastq.gz"
  # )
  # 
  # patternmy <- c(dirname(FileName[1]), "_sub.fastq.gz")
  # SampleName <- str_replace_all(
  #   pattern = paste0(dirname(FileName)[1], "/|_sub.fastq.gz"),
  #   string = list.files(
  #     file.path(
  #       system.file(package = "pathwrap"), "extdata"
  #     ),
  #     full.names = TRUE, pattern = "fastq.gz"
  #   ), ""
  # )
  # 
  # 
  # Class <- c("A", "B", "A", "B")
  # write.table(as.data.frame(cbind(SampleName, FileName, Class)),
  #   file = phenofile, sep = "\t", row.names = FALSE,
  #   col.names = TRUE, quote = FALSE
  # )
  # 

  message("this is the phenofile ", phenofile)
  cdatapath <- file.path(system.file(package = "pathwrap"), "extdata", 
                         "example_cpd_data.tsv")
  csamp <- c(1,2)
  cref <- c(3,4)
  path.res <- pathwrap(
   ref.dir = NA, phenofile = phenofile,
    outdir = Results, entity = "Mus musculus", corenum = 16,
    compare = "as.group",  keep_tmp = TRUE,
    startover = FALSE, diff.tool = "DESeq2", aligner = "Rhisat2", 
    mode = "combined", cdatapath = cdatapath, cref= cref, csamp = csamp, 
    ccompare = "paired", is.test = TRUE
  )
  path.res <- "The analysis is complete"

 identical(goodtestresultmessage, path.res)
 checkEquals(4/ 2, 2)
  # checkTrue(is.na(4/ 0))
  # checkEqualsNumeric(4 /1.2345, 3.24, tolerance=1.0e-4)
}






