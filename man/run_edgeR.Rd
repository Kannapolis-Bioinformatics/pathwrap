% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_rna_edgeR.R
\name{run_edgeR}
\alias{run_edgeR}
\title{Run standard edgeR2 for differential gene expression analysis and
plot volcano plots}
\usage{
run_edgeR(cnts, grp.idx, edger.dir)
}
\arguments{
\item{cnts}{: counts of genes}

\item{grp.idx}{: index of the reference and sample for differential analysis}

\item{edger.dir}{: the directory in which edgeR results will be stored}
}
\value{
log fold expression values
}
\description{
run_edgeR takes counts and the list indicating reference and samples and
the directory where the results are stored and performs the edgeR analysis
The output is result table with columns of genes and logFC from result of
edgeR analysis and a volcanoplot.
It returns the logFC to be used by GAGE for gene set analysis.
}
