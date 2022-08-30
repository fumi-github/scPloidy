#' Liver Cells from Rat
#'
#' The dataset includes 3572 nuclei from the liver of a male rat of 16 weeks age,
#' which was fed normal diet.
#' Overlapping of single-nucleus ATAC-seq fragments was computed with the
#' `fragmentoverlapcount` function and saved as `fragmentoverlap`.
#' The cell type of the nuclei are saved in the dataframe `cells`.
#' The data for rat SHR_m154211 was taken from publication.
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
