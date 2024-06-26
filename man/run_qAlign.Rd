% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_rna_qalign.R
\name{run_qAlign}
\alias{run_qAlign}
\title{\enumerate{
\item RUN THE ANALYSIS # Alignment and counting
}}
\usage{
run_qAlign(
  corenum,
  endness,
  sampleFile,
  genomeFile,
  geneAnnotation,
  ref.dir,
  cacheDir,
  aligner
)
}
\arguments{
\item{corenum}{: the number of cores used during alignment}

\item{endness}{: weather its paired end or single end}

\item{sampleFile}{: the file where information about location of sample is
stored, see qAlign for more}

\item{genomeFile}{: the genome file used for alignment or bioconductor genome
package}

\item{geneAnnotation}{: gene annnotaion file used for gene counting or
bioconductor annotation package}

\item{ref.dir}{: directory in which genomeFile and genomeAnnotaion are stored}

\item{cacheDir}{: directory where temporary files generated during alignment
are store}

\item{aligner}{: weather Rhisat2 or Rbowtie should be used for alignment}
}
\value{
R object generated from the alignment step
}
\description{
this function runs the Rhisat2 or Rbowtie for alignment on paired or
single end mode. It saves the alignment object in RDS file which can
be loaded in R for further analysis. If the reference index is not found
in the reference directory, it creates reference index before running
alignment. If the reference genome is a package the reference index is
created as R package.It generates the barplot of mapped and unmapped sequence
reads.
}
