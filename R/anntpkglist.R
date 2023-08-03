#' Genome and Annotaion packages in Bioconductor
#'
#' Chart for latest genome and annotation package for few selected species that
#' are downloaded when no reference is provided by user
#' @return the data frame with packages to install before analysis
#' @format ## `anntpkglist` A dataframe with 32 rows and three columns;
#' rows describe different species
#' \describe{
#'   \item{genome}{Name of the genome package}
#'   \item{species}{Scientific names of species}
#'   \item{annotation}{Name of the annotation species}
#'   ...
#' }
#' @source <https://www.bioconductor.org/packages/release/data/annotation/>
"anntpkglist"
