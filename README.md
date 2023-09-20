# pathwrap, wrapper for pathview

## Overview

Pathwrap is a RNASeq analysis tool that provides a wrapper for the processing of RNAseq datasets from quality control of raw reads to the visualization of enriched pathwys obtained from processing the datasets. This tool runs all the essential steps of RNAseq processing from quality control, filtering out low quality reads, trimming adapters, sequence alignemnt, alignment count, differential analysis, enrichemnt analysis and pathway visualization. It is the first tool that combines all essential steps of RNSeq analysis till pathway visualization. It has the ability to continue the analysis if it is halted at any stage and generate quality pictures and generate comprehensive analysis of the data. 

## Installation
In order to install pathwrap, open R (version "4.3") and write

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Kannapolis-Bioinformatics/pathwrap", force = TRUE, build_vignette = TRUE)
```

Also you can find the latest annotation and genome package useful for analysis by running following code.

```r 
library(pathviewwrap)
data(anntpkglist)
genomepkg <- anntpkglist$genome[which(anntpkglist$species=="Homo sapiens")]
anntpkg <- anntpkglist$annotation[which(anntpkglist$species=="Homo sapiens")]
BiocManager::install(genomepkg)
BiocManager::install(anntpkg)
```

## Quick start with demo data 
Just run the pathwrap function with as much argument as possible for compelte analysis. You will need phenofile which has information about the path in which the raw files are stored and the class or category each sample belong to.

``` r
# This code creates the phenofile and runs the pathviewrap
#this is a demo and phenofile can be created in any way.

#create directory to store results
Results <- tempdir()
#Make sure research is path to the location where you can see
#the data and explore it ; like
#>Results <- "/Users/edhungel/Research/Documents/myresults
  
#phenofile should be path to some file not temporary file
phenofile <-tempfile("hellotmpphenofile.txt")
#Make sure this is a file path readable by R, read.table like 
#>phenofile <- "/Users/edhungel/Research/Documents/myphenofile.txt"

#create columns for phenofile, this is for SE data
#col.names should be SampleName, FileName and Class for SE data
library(stringr)
FileName <- list.files(file.path(system.file(
    package = "pathviewwrap"), "extdata"), pattern = "fastq.gz",
    full.names = TRUE)

#patternmy <- c(dirname( FileName[1]) , "_sub.fastq.gz")
SampleName <- str_replace_all(
    pattern = paste0(dirname(FileName)[1], "/|_sub.fastq.gz" ),
    string = list.files(file.path(
    system.file(package = "pathviewwrap"), "extdata"), 
    full.names = TRUE,
    pattern = "fastq.gz" ) , "" )


Class <- c("A", "B", "A", "B")
write.table(as.data.frame(cbind(SampleName, FileName, Class)), 
            file = phenofile, sep = "\t", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)


message("this is the phenofile ", phenofile )
library(pathviewwrap)
system.time({
    pathviewwrap(
        ref.dir = NA, phenofile = phenofile,
        outdir = Results, entity = "Mus musculus", corenum = 16,
        compare = "as.group", seq_tech = "Illumina", keep_tmp = TRUE,
        startover = TRUE, diff.tool = "DESeq2", aligner = "Rhisat2"
    )
})

```

## Steps run by the wrapper 
The steps run by the wrapper are as follows:

#test

With one wrapper function, it runs all the steps listed below. 

## STEP 1 : Quality control

# STEP 1a: running fastqc

It runs fastqc analysis in R using fastqcr. If fastqc is not available in system to run by R, this function is capable of downloading the fastqc tools before running the quality check. The results are standard html files where the quality of each fastq files can be examined.

# STEP 1b : running fastp

After running fastqc, it runs fastp. The function takes name of the samples and for each sample does the quality and adapter trimming for Illumina and long read sequencing. It works for both PE and SE data. HTML files are generated for each fastq files that has information/figures of quality control before and after quality trimming.

## STEP 2: making txdb obj

The wrapper then makes TxDb object either from annnotation file or by loading from the annotation package. Make_txdbobj uses the makeTxDbFromGFF to make TxDb object from transcript annotations available in gtf file in ref.dir or if ref.dir is NA, it makes txdb object fromt the annotaion package

## STEP 3 : Alignmnet and counting 2 

After the txdb object is formed, the wrapper runs the Rhisat2 or Rbowtie for alignment on paired or single end mode depending on data. It saves the alignment object in RDS file which can be loaded in R for further analysis. If the reference index is not found in the reference directory, it creates reference index before running alignment. If the refernece genome is a package the reference index is created as R package. It generates the barplot of mapped and unmapped sequence reads.

## STEP 4: counting aligned sequences

After aligning the reads to reference genome, the wrapper generates the count of gene and store it in a table with genes in row and counts in columns. The table is stored as a RDS file that can be loaded into R for future analysis and if gene id is ensembl, the wrapper converts it to entrez to match with genes names for gene set analysis. 

# Differential gene analysis

## STEP 5a ; running differential gene analysis using DESeq2

Then the wrapper runs standard DESeq2 for differential gene expression analysis and plots volcano plots. The function run_deseq2 takes counts and the list indicating reference and samples and the directory where the results are stored and performs the deseq2 analysis. The output is result table with columns of genes and log2FoldChange from result of deseq2 analysis and a volcanoplot.


## STEP 5b ; running differential gene analysis using edgeR

Then the wrapper runs standard edgeR for differential gene expression analysis and plots volcano plots. The function run_deseq2 takes counts and the list indicating reference and samples and the directory where the results are stored and performs the deseq2 analysis. The output is result table with columns of genes and log2FoldChange from result of deseq2 analysis and a volcanoplot.

## STEP 6 : running pathway analysis using GAGE 

After the differential gene analysis the wrapper runs generally applicable gene set enrichment for pathway analysis, GAGE based upon the user supplied comparision method for the species specified. The biological process, cellular component and molecular function analysis for GO terms are done seperately. Also, KEGG disease and KEGG signalling and metabolism pathways are analysed seperately. 

## STEP 7: visualizing the pathway using Pathview

Finally the top enriched pathways with "q.val" < 0.1 are visualized using pathview.

## More information
Please watch out for paper in making. 
https://docs.google.com/document/d/1pfMI-umnS7GCW9aoAqEm0tZv9g6eVKSA/edit


Thank you for your interest.

Please send all queries to [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)<br />
[Eliza Dhungel](mailto:edhungel@uncc.ed) <br /> 
Or [open an issue]
(https://github.com/Kannapolis-Bioinformatics/pathwrap/issues)

