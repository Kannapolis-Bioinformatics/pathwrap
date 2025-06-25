# Pathwrap

## Overview

Pathwrap is an analysis tool for the processing of RNAseq datasets from raw data to data visualizations. Pathwrap is built on pathway enrichment tool GAGE (Generally Applicable Gene-set Enrichment for Pathway Analysis) and pathway visualization using Pathview.  Features include all the essential steps of RNAseq processing including read quality control (e.g., trimming and filtering), read mapping,  read summarization/quantification, statistical differential abundance analysis (DESeq2 and edgeR), pathway enrichment (GAGE using KEGG KO), and pathway visualization (Pathview). Pathwrap provides a start to finish automatic pipeline within the R framework for comprehensive analysis of RNAseq data. In addition it allows seamless integration of pathway analysis and visualization of RNAseq data with quantitative metabolomics data.

## Installation
1. In order to install pathwrap, open R (version "4.4") and write

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("remotes")
BiocManager::install("Kannapolis-Bioinformatics/pathwrap@simplified_beta")
```

2. You can find the latest annotation and genome package useful for analysis by running following code.

```r 
library(pathwrap)
data(anntpkglist)
#change Mus musculus to Homo sapiens if you are analysing human reads
genomepkg <- anntpkglist$genome[which(anntpkglist$species=="Mus musculus")]
anntpkg <- anntpkglist$annotation[which(anntpkglist$species=="Mus musculus")]

#run codes below to install the packages
#this is necessary if path to reference directory is not provided
#BiocManager::install(genomepkg)
#BiocManager::install(anntpkg)
```

If path to reference directory is provided using argument ref.dir, the path should be writable to create indexes and it should contain both genome reference (.fa,/.fasta) and genome annotaion file (.gtf/.gff). 

3. A phenofile is necessary for any run of pathwrap. The phenofile should be tab delimited file with information about the path of raw files and class to which each sample belong to. Phenofile for single ended reads looks like this.
```{r sampleFileSingle, echo=FALSE, results='asis'}
SampleName	FileName	Class
sample1	/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library/pathwrap/extdata/sample1_sub.fastq.gz	A
sample3	/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library/pathwrap/extdata/sample3_sub.fastq.gz	B
sample5	/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library/pathwrap/extdata/sample5_sub.fastq.gz	A
sample6	/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library/pathwrap/extdata/sample6_sub.fastq.gz	B

```
Phenofile for paired end reads have at least 4 columns with column names as: SampleName    FileName1    FileName2	Class.

In case of paired experiment design, the phenofile should have extra column named PairedInfo to indicate sample pairs. 


## Quick start with demo data 
To run pathwrap, minimum required arguments are path to phenofile and scientic name of the species of interest. 

```r
library(pathwrap)
pathwrap(phenofile=file.path(system.file(package = "pathwrap"), "extdata", "phenofile_SE.txt"),
                    entity= "Mus musculus")
```

An example of how phenofile can be created is as follows. 
``` r
# This code creates the phenofile and runs the wrapper for Pathview
#this is a demo and phenofile can be created in any way.

#create directory to store results
Results <- tempdir()
#Make sure results is path to the location where you can see
#the data and explore it ; like
#>Results <- "/Users/edhungel/Research/Documents/myresults"
  
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
#is the sample reference or experiment? make sure this matches the row 
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
csamp <- c(1,2)
cref <- c(3,4)
if(interactive()){ system.time({
    pathwrap(
     phenofile = phenofile,
    entity = "Mus musculus", 
    cdatapath=cdatapath, 
    csamp=csamp,cref=cref
    )
})}

```
## Parameter to run pathwrap 
## Required
| Parameter | Description | Example |
|--------|--------|--------|
|phenofile | file where the path of raw data is stored | "/usr/document/myrawfile" |
|entity | Scientific name of species whose RNA is being analyzed | "Homo sapiens " |

## Optional

| Parameter | Description | Default  | 
|--------|--------|--------|
| outdir | main directory for storing output of the process | "./Results" |
| startover |  do you want to start from beginning, | FALSE |
|corenum  | number of cores avaialble for run | detectCores() |
| ref.dir |path to reference directory which contain reference file(*.fa) and annotation file(*.gtf) | NA |
| cacheDir | directory where temporary files created during alignment | tempdir() |
|aligner |One of "Rhisat2" or "Rbowtie2"; Rbowtie2 can be very slow for human and eukaryotic species | "Rhisat2" |
| gcompare |how the comparision is done for transcripts/genes | "unpaired"  | 
| npca |number of genes to use for pca | 19 |
|nheatmap |number of genes for heatmap | 10 |
|cdatapath |data path for compound data | NA |
|cpd_id_type | "KEGG COMPOUND accession" | "KEGG COMPOUND accession" |
|csamp |index/row number where sample files are ex: c(5,6,4) | NULL |
|cref |index/row number where references are, ex: c(1,2,3) | NULL  |
|ccompare | how the compound data is compared |   "paired" | 
|qcut |threshold for pathway selection | 0.01 |
|pathids |pathway of interest only necessary if enrichment is not run | "04110" |
|nchunks | default 1, in how many chunks you want to run alignment | 1 |
|keep_tmp |weather to store aligned bam files and trimmed fastq files | TRUE | 
|diff.tool | weather to use "DESeq2" or edgeR for differential gene analysis | "DESeq2" |



## Steps run by the pathwrap
The steps run are as follows:

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

Or

## STEP 5b ; running differential gene analysis using edgeR

Then the wrapper runs standard edgeR for differential gene expression analysis and plots volcano plots. The function run_deseq2 takes counts and the list indicating reference and samples and the directory where the results are stored and performs the DESeq2 analysis. The output is result table with columns of genes and log2FoldChange from result of DESeq2 analysis and a volcanoplot.

## STEP 6 : running pathway analysis using GAGE 

After the differential gene analysis the wrapper runs generally applicable gene set enrichment for pathway analysis, GAGE based upon the user supplied comparision method for the species specified. The biological process, cellular component and molecular function analysis for GO terms are done seperately. Also, KEGG disease and KEGG signalling and metabolism pathways are analysed seperately. 

## STEP 7: visualizing the pathway using Pathview

Finally the top enriched pathways with "q.val" < 0.01 are visualized using Pathview.

## More information
Please watch out for paper in making. 
https://docs.google.com/document/d/1pfMI-umnS7GCW9aoAqEm0tZv9g6eVKSA/edit


Thank you for your interest.

Please send all queries to [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)<br />
[Eliza Dhungel](mailto:edhungel@uncc.ed) <br /> 
Or [open an issue]
(https://github.com/Kannapolis-Bioinformatics/pathwrap/issues)

