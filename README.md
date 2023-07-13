# pathwrap tool
Pathwrap is a RNASeq analysis tool that provides a wrapper for the processing of RNAseq datasets from quality control of raw reads to the visualization of enriched pathwys obtained from processing the datasets. This tool runs all the essential steps of RNAseq processing from quality control, filtering out low quality reads, trimming adapters, sequence alignemnt, alignment count, differential analysis, enrichemnt analysis and pathway visualization. It is the first tool that combines all essential steps of RNSeq analysis till pathway visualization. 

In order to install pathwrap
```
# Pre-installed BiocManager
BiocManager::install("Kannapolis-Bioinformatics/pathwrap", ref="rnaseqwrap", force = T)
``

Please send all queries to [Eliza Dhungel](mailto:edhungel@uncc.ed) <br />  
[Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)<br />
Or [open an issue](https://github.com/Kannapolis-Bioinformatics/pathwrap/issues)
