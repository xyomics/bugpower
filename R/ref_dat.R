#' Reference Data
#'
#' The data set contains the mean distributions and standard deviation quartiles of
#' species abundance, metagenomic pathway abundance and metatranscriptomic pathway abundance
#' (normalized by metagenomics data).
#'
#' @docType data
#'
#' @usage data(ref_data)
#'
#' @format A list contains mean distribution vectors and quartile matrix for each type of data.
#' \describe{
#'   \item{mean_distribution_bugs}{A vector of mean distribution of species relative abundance.}
#'   \item{sd_quartile_bugs}{A vector of 3 quartiles of species relative abundance.}
#'   \item{mean_distribution_pathway_dna}{A vector of relative abundance of metagenomic pathways.}
#'   \item{sd_quartile_pathway_dna}{A vector of 3 quartiles of relative abundance of metagenomic pathways.}
#'   \item{mean_distribution_pathway_rna}{A vector of relative abundance of metatranscritomic pathways.}
#'   \item{sd_quartile_pathway_rna}{A vector of 3 quartiles of relative abundance of metatranscriptomic pathways.}
#' }
#'
#' @keywords datasets
#'
#' @references STAR
#'
#' @source \href{http://xyomics}{STAR data arcive}
#'
#' @examples
#' data(ref_dat)
"ref_dat"

