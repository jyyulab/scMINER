#' Raw count matrix of PBMC14k dataset
#'
#' This dataset contains the raw UMI counts of 14,000 Peripheral Blood Mononuclear Cells (PBMCs) of 7 known cell types, with 2,000 cell each.
#'
#' @format ## `pbmc14k_rawCount`
#' A large dgCMatrix with 17,986 rows and 14,000 columns:
#' \describe{
#'   \item{CACTTTGACGCAAT}{UMI counts of all genes of the cell "CACTTTGACGCAAT"}
#'   \item{GTTACGGAAACGAA}{UMI counts of all genes of the cell "GTTACGGAAACGAA"}
#'   \item{AGTCACGACAGGAG}{UMI counts of all genes of the cell "AGTCACGACAGGAG"}
#'   ...
#' }
#' @source A subset of Filtered_DownSampled_SortedPBMC_data.csv from <https://zenodo.org/record/3357167#.YhQNF2RKj6V>
"pbmc14k_rawCount"
