#' Liver Cells from a Rat
#'
#' The dataset includes 3572 nuclei obtained from the liver of
#' a 16 weeks old male rat, which was fed normal diet.
#' Overlapping of single-nucleus ATAC-seq fragments was computed with the
#' \code{fragmentoverlapcount} function and saved as \code{fragmentoverlap}.
#' The cell type of the nuclei are saved in the data.frame \code{cells}.
#' The data for rat SHR_m154211 was taken from the publication cited below.
#'
#' @docType data
#'
#' @usage data(SHR_m154211)
#'
#' @source \href{https://doi.org/10.1101/2022.07.12.499681}{Takeuchi et al. (2022) bioRxiv}
#'
#' @examples
#' data(SHR_m154211)
#' fragmentoverlap = SHR_m154211$fragmentoverlap
#' p = ploidy(fragmentoverlap, c(2, 4, 8))
#' head(p)
#' cells = SHR_m154211$cells
#' table(cells$celltype, p$ploidy.moment[match(cells$barcode, p$barcode)])
"SHR_m154211"
