#' Species-by-site raw incidence data
#'
#' This data set consists of the species-by-site raw incidence (detection/non-detection) matrix and the
#' phylogenetic tree of 41 bird species collected in November 2012 at Barrington Tops National Park,
#' Australia from two sites: the North Site (27 species from \code{nT =} 12 point counts) and the South
#' Site (38 species from \code{nT =} 17 point counts). Each point count is regarded as a sampling unit in
#' which species incidence (detection or non-detection) is recorded. See the vignette for part of the
#' data in the required input format, and Chao et al. (2015) for the statistical analysis of this data set.
#'
#' @usage data(data.inc)
#' @format A list includes three objects:
#' \describe{
#'   \item{$data}{a data frame with species-by-site raw incidence records for 41 species.}
#'   \item{$tree}{a phylo object giving the phylogenetic tree in Newick format.}
#'   \item{$nT}{a vector giving the numbers of sampling units in all sites.}
#' }
#' @source Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015). Rarefaction and extrapolation of phylogenetic diversity. \emph{Methods in Ecology and Evolution}, 6, 380-388.
"data.inc"
