utils::globalVariables(c("FileName", "filenames", "endness", "trim.dir", "seq_tech",
                         "sig_n_met","edgeR.dir", "disease", "signalling", "metabolism", "korg","bods",
                         "biological_process", "molecular_function","cellular_component","aligned_bam"))



utils::globalVariables(c("fastp_log","tot.seq", "module", "status", "fastqc_results", "differential_analysis", "fastp_results", "gage_results", "pathway_analysis","edgeR","DESeq2", "KEGG", "GO" , "anntpkglist"))


# DESeq2 DataFrame FileName GO KEGG aligned_bam anntpkglist
# biological_process bods cellular_component differential_analysis
# disease edgeR endness fastp_log fastp_results fastqc_results
# filenames gage_results korg metabolism module molecular_function
# pathway_analysis renameSeqlevels seq_tech seqlevels seqnames
# sig_n_met signalling status str_remove tot.seq trim.dir width
