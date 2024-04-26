#' Basal cell carcinoma sample SU008_Tumor_Pre
#'
#' The dataset includes 788 nuclei obtained from
#' basal cell carcinoma sample SU008_Tumor_Pre.
#' Overlapping of single-nucleus ATAC-seq fragments was computed with the
#' \code{fragmentoverlapcount} function.
#'
#' @name GSE129785_SU008_Tumor_Pre
#' @docType data
#' @usage data(GSE129785_SU008_Tumor_Pre)
#' @references Satpathy et al. (2019) Nature Biotechnology 37:925 \doi{10.1038/s41587-019-0206-z}
#' @source \href{https://www.ncbi.nlm.nih.gov/geo}{GEO, GSE129785}
#'
#' @examples
#' data(GSE129785_SU008_Tumor_Pre)
#' levels = c(2, 4)
#' result = cnv(SU008_Tumor_Pre_fragmentoverlap,
#'              SU008_Tumor_Pre_windowcovariates,
#'              levels = levels,
#'              deltaBICthreshold = -600)

#' @rdname GSE129785_SU008_Tumor_Pre
#' @format \code{SU008_Tumor_Pre_fragmentoverlap} is a dataframe of fragmentoverlap.
"SU008_Tumor_Pre_fragmentoverlap"

#' @rdname GSE129785_SU008_Tumor_Pre
#' @format \code{SU008_Tumor_Pre_windowcovariates} is a dataframe of windows and peaks.
"SU008_Tumor_Pre_windowcovariates"

#' @rdname GSE129785_SU008_Tumor_Pre
#' @format \code{rescnv} is a list containing the output of cnv function.
"rescnv"
