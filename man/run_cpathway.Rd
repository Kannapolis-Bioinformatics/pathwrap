% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_cpathway_analysis.R
\name{run_cpathway}
\alias{run_cpathway}
\title{Function for running the compound pathway analysis}
\usage{
run_cpathway(
  cdatapath,
  cpd_id_type = "KEGG COMPOUND accession",
  csamp,
  cref,
  ccompare = "paired",
  cset_dir,
  entity
)
}
\arguments{
\item{cdatapath}{:pathway for compound data}

\item{cpd_id_type}{type of compound id as compatible to ones in
pathview::data(rn.list)}

\item{csamp}{column number of sample}

\item{cref}{column number of reference}

\item{ccompare}{comparision type for sample}

\item{cset_dir}{directory where results of compound set analysis is stored}

\item{entity}{scientific name of the organism of study}
}
\value{
gage_return_obj_c return object from run_gage for compound
}
\description{
This function prepares for running the gage analysis on compound sets
}
