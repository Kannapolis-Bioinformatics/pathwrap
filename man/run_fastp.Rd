% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_rna_fastp.R
\name{run_fastp}
\alias{run_fastp}
\title{Runs fastp for quality and adapter trimming}
\usage{
run_fastp(sampleName, FileName, endness, trim.dir, corenum)
}
\arguments{
\item{sampleName}{: name of the sampleverbiage}

\item{FileName}{name of the raw files, dataframe with colnames
FileName1 and FileName2 for paired read}

\item{endness}{wheatehr the raw data is paired or single ended}

\item{trim.dir}{directory to store the trimmed fasta files}

\item{corenum}{number of threads}
}
\value{
no value returned
}
\description{
This function takes name of the samples and for each sample does the quality
and adapter trimming for Illumina and long read sequencing. It works for both
PE and SE data
}
