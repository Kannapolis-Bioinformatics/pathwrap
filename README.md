# Pathwrap

## Overview

Pathwrap is an analysis tool for the processing of RNAseq datasets from raw data to data visualizations. Pathwrap is built on pathway enrichment tool GAGE (Generally Applicable Gene-set Enrichment for Pathway Analysis) and pathway visualization using pathview.  Features include all the essential steps of RNAseq processing including read quality control (e.g., trimming and filtering), read mapping,  read summarization/quantification, statistical differential abundance analysis (DESeq2 and edgeR), pathway enrichment (GAGE using KEGG KO), and pathway visualization (pathview). Pathwrap provides a start to finish automatic pipeline within the R framework for comprehensive analysis of RNAseq data. In addition it allows seamless integration of pathway analysis and visualization of RNAseq data with quantitative metabolomics data.

## Installation
In order to install pathwrap, open R (version "4.3") and write

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("pathwrap")
```

Also you can find the latest annotation and genome package useful for analysis by running following code.

```r 
library(pathwrap)
data(anntpkglist)
genomepkg <- anntpkglist$genome[which(anntpkglist$species=="Mus musculus")]
anntpkg <- anntpkglist$annotation[which(anntpkglist$species=="Mus musculus")]
BiocManager::install(genomepkg)
BiocManager::install(anntpkg)
```

## Quick start with demo data 
Just run the pathwrap function with as much argument as possible for complete analysis. You will need phenofile which has information about the path in which the raw files are stored and the class or category each sample belong to.

``` r
# This code creates the phenofile and runs the wrapper for pathview
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
    package = "pathwrap"), "extdata"), pattern = "fastq.gz",
    full.names = TRUE)

SampleName <-str_remove_all( basename(FileName), ".fastq.gz")
#patternmy <- c(dirname( FileName[1]) , "_sub.fastq.gz")

Class <- c("A", "B", "A", "B")
write.table(as.data.frame(cbind(SampleName, FileName, Class)), 
            file = phenofile, sep = "\t", row.names = FALSE, 
            col.names = TRUE, quote = FALSE)
cdatapath <- file.path(system.file(package = "pathwrap"), "extdata", 
                         "example_cpd_data.tsv")

message("this is the phenofile ", phenofile )
library(pathwrap)
system.time({
    pathwrap(
            ref.dir = NA, phenofile = phenofile,mode = "gene", 
    outdir = Results, entity = "Mus musculus", corenum = 16,
    compare = "as.group",  keep_tmp = TRUE,
    startover = TRUE, diff.tool = "DESeq2", aligner = "Rhisat2",
    cdatapath = cpath, cref= cref, csamp = csamp, ccompare = "paired"
    )
})

```

## Steps run by the pathwrap
The steps run are as follows:

#test

With one function, it runs all the steps listed below. 

## STEP 1 : Quality control

# STEP 1a: running fastqc

It runs fastqc analysis in R using fastqcR. If fastqc is not available in system to run by R, this function is capable of downloading the fastqc tools before running the quality check. The results are standard html files where the quality of each fastq files can be examined.

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

Then the wrapper runs standard edgeR for differential gene expression analysis and plots volcano plots. The function run_deseq2 takes counts and the list indicating reference and samples and the directory where the results are stored and performs the DESeq2 analysis. The output is result table with columns of genes and log2FoldChange from result of DESeq2 analysis and a volcanoplot.

## STEP 6 : running pathway analysis using GAGE 

After the differential gene analysis the wrapper runs generally applicable gene set enrichment for pathway analysis, GAGE based upon the user supplied comparision method for the species specified. The biological process, cellular component and molecular function analysis for GO terms are done seperately. Also, KEGG disease and KEGG signalling and metabolism pathways are analysed seperately. 

## STEP 7: visualizing the pathway using Pathview

Finally the top enriched pathways with "q.val" < 0.01 are visualized using pathview.

## More information
Please watch out for paper in making. 
https://docs.google.com/document/d/1pfMI-umnS7GCW9aoAqEm0tZv9g6eVKSA/edit


Thank you for your interest.

Please send all queries to [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)<br />
[Eliza Dhungel](mailto:edhungel@uncc.ed) <br /> 
Or [open an issue]
(https://github.com/Kannapolis-Bioinformatics/pathwrap/issues)

