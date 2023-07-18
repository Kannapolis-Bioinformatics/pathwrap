# pathwrap tool
Pathwrap is a RNASeq analysis tool that provides a wrapper for the processing of RNAseq datasets from quality control of raw reads to the visualization of enriched pathwys obtained from processing the datasets. This tool runs all the essential steps of RNAseq processing from quality control, filtering out low quality reads, trimming adapters, sequence alignemnt, alignment count, differential analysis, enrichemnt analysis and pathway visualization. It is the first tool that combines all essential steps of RNSeq analysis till pathway visualization. 

In order to install pathwrap
BiocManager::install("raw-lab/pathwrap", ref="rnaseqwrap", force = T) #make sure BiocManager is already installed within the R 

#test
