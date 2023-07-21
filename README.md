# pathwrap tool
Pathwrap is a RNASeq analysis tool that provides a wrapper for the processing of RNAseq datasets from quality control of raw reads to the visualization of enriched pathwys obtained from processing the datasets. This tool runs all the essential steps of RNAseq processing from quality control, filtering out low quality reads, trimming adapters, sequence alignemnt, alignment count, differential analysis, enrichemnt analysis and pathway visualization. It is the first tool that combines all essential steps of RNSeq analysis till pathway visualization. 

It has the ability to continue the analysis if it is halted at any stage and generate quality pictures and generate comprehensive analysis of the data. 

In order to install pathwrap
BiocManager::install("raw-lab/pathwrap", ref="rnaseqwrap", force = T) #make sure BiocManager is already installed within the R 

#test

With one wrapper function, it runs all the steps listed below. 

STEP 1 : Quality control

STEP 1a: running fastqc

It runs fastqc analysis in R using fastqcr. If fastqc is not available in system to run by R, this function is capable of downloading the fastqc tools before running the quality check. The results are standard html files where the quality of each fastq files can be examined.

STEP 1b : running fastp

After running fastqc, it runs fastp. The function takes name of the samples and for each sample does the quality and adapter trimming for Illumina and long read sequencing. It works for both PE and SE data. HTML files are generated for each fastq files that has information/figures of quality control before and after quality trimming.

STEP 2: making txdb obj

The wrapper then makes TxDb object either from annnotation file or by loading from the annotation package. Make_txdbobj uses the makeTxDbFromGFF to make TxDb object from transcript annotations available in gtf file in ref.dir or if ref.dir is NA, it makes txdb object fromt the annotaion package

STEP 3 : Alignmnet and counting 2 

After the txdb object is formed, the wrapper runs the Rhisat2 or Rbowtie for alignment on paired or single end mode depending on data. It saves the alignment object in RDS file which can be loaded in R for further analysis. If the reference index is not found in the reference directory, it creates reference index before running alignment. If the refernece genome is a package the reference index is created as R package. It generates the barplot of mapped and unmapped sequence reads.

STEP 4: counting aligned sequences

After aligning the reads to reference genome, the wrapper generates the count of gene and store it in a table with genes in row and counts in columns. The table is stored as a RDS file that can be loaded into R for future analysis and if gene id is ensembl, the wrapper converts it to entrez to match with genes names for gene set analysis. 

STEP 5a ; running differential gene analysis using DESeq2

Then the wrapper runs standard DESeq2 for differential gene expression analysis and plots volcano plots. The function run_deseq2 takes counts and the list indicating reference and samples and the directory where the results are stored and performs the deseq2 analysis. The output is result table with columns of genes and log2FoldChange from result of deseq2 analysis and a volcanoplot.


STEP 5b ; running differential gene analysis using edgeR

Then the wrapper runs standard edgeR for differential gene expression analysis and plots volcano plots. The function run_deseq2 takes counts and the list indicating reference and samples and the directory where the results are stored and performs the deseq2 analysis. The output is result table with columns of genes and log2FoldChange from result of deseq2 analysis and a volcanoplot.

STEP 6 : running pathway analysis using GAGE 

After the differential gene analysis the wrapper runs generally applicable gene set enrichment for pathway analysis, GAGE based upon the user supplied comparision method for the species specified. The biological process, cellular component and molecular function analysis for GO terms are done seperately. Also, KEGG disease and KEGG signalling and metabolism pathways are analysed seperately. 

STEP 7: visualizing the pathway using Pathview

Finally the top enriched pathways with  "q.val" < 0.1 are visualized using pathview.

