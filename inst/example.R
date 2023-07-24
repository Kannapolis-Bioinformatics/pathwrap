#here write the code to run for the test of mouse , use this file and the 
#phenofile in the package
#BiocManager::install("Kannapolis-Bioinformatics/pathwrap", force  = T)
Results <-tempdir() 
library(pathviewwrap)
system.time({
  pathviewwrap(ref.dir = NA , 
               phenofile= system.file(package= "pathviewwrap",
                                      "extdata/phenotest_ungz.txt") ,
               outdir=Results, entity="Mus musculus", corenum = 16, 
               compare = "as.group", seq_tech="Illumina",keep_tmp = TRUE, 
               rerun = FALSE, diff.tool = "DESeq2", aligner = "Rhisat2")
})


