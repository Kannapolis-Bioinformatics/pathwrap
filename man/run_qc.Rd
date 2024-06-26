% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_rna_fastqc.R
\name{run_qc}
\alias{run_qc}
\title{Run fastqc analysis in R}
\usage{
run_qc(fq.dir, qc.dir, corenum)
}
\arguments{
\item{fq.dir}{: the directory in which raw RNAseq files are stored}

\item{qc.dir}{: the directory in which results of quality check are stored}

\item{corenum}{: the number of cores used for running quality check}
}
\value{
message of fastqc
}
\description{
runs fastqc analysis in R using fastqcr. If fastqc is not available in system
to run by R. this function is capable of downloading the fastqc tools before
running the quality check. The results of quality check is aggregated and
barplot of the total sequence and heatmap of the status of the qc check is
produced. do not install fastqc in your results directory
}
