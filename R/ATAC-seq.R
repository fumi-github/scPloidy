#' Count Overlap of ATAC-seq Fragments
#'
#' Count Overlap of ATAC-seq Fragments
#'
#' @param file Filename of the file for ATAC-seq fragments.
#' The file must be block gzipped (using the \code{bgzip} command)
#' and accompanied with the index file (made using the \code{tabix} command).
#' The uncompressed file must be a tab delimited file,
#' where each row represents one fragment.
#' The four columns are chromosome name, start position, end position,
#' and barcode (i.e., name) of the cell including the fragment.
#' See vignette for details.
#' @param targetregions GRanges object for the regions where overlaps are counted.
#' Usually all of the autosomes.
#' @param excluderegions GRanges object for the regions to be excluded.
#' Simple repeats in the genome should be listed here,
#' because repeats can cause false overlaps.
#' A fragment is discarded if its 5' or 3' is located in \code{excluderegions}.
#' @param targetbarcodes Character vector for the barcodes of cells to be analyzed,
#' such as those passing quality control.
#' If \code{NULL}, all barcodes in the input file are analyzed.
#' @return A tibble with each row corresponding to a cell.
#' For each cell, its barcode, the total count of the fragments \code{nfrag},
#' and the count distinguished by overlap depth are given.
#'
#' @importFrom dplyr arrange desc distinct group_by mutate n rename summarize ungroup
#' @importFrom GenomicRanges findOverlaps makeGRangesFromDataFrame
#' @importFrom Rsamtools scanTabix TabixFile
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar
#' @export
fragmentoverlapcount = function (file,
                                 targetregions,
                                 excluderegions = NULL,
                                 targetbarcodes = NULL) {
  sumoverlaplist = list()
  tbx = TabixFile(file = file)
  pb = txtProgressBar(max = length(targetregions), style = 3)
  for (i in 1:length(targetregions)) {

    # Load fragments.
    res = scanTabix(tbx, param = targetregions[i])
    frags = read.csv(textConnection(res[[1]]),
                     sep = "\t",
                     header = FALSE)
    colnames(frags) = c("chr", "start", "end", "BC")

    if (! is.null(targetbarcodes)) {
      frags = frags[frags$BC %in% targetbarcodes, ]
    }

    # Discard "semi-duplicate" fragments;
    # Hypothesized Tn5 transposition only at one strand.
    frags = frags %>%
      group_by(BC) %>%
      # Between chrX:100-200 and chrX:100-300, only retain first
      arrange(start, end) %>%
      distinct(start, .keep_all = TRUE) %>%
      # Between chrX:100-300 and chrX:200-300, only retain last
      arrange(desc(end), desc(start)) %>%
      distinct(end, .keep_all = TRUE) %>%
      ungroup() %>%
      arrange(start, end, BC)

    # Discard fragment if 5' or 3' is located in excluderegions.
    if (! is.null(excluderegions)) {
      query = makeGRangesFromDataFrame(
        data.frame(
          seqnames = frags$chr,
          start    = frags$start,
          end      = frags$start))
      x = is.na(findOverlaps(query, excluderegions, select = "first"))
      query = makeGRangesFromDataFrame(
        data.frame(
          seqnames = frags$chr,
          start    = frags$end,
          end      = frags$end))
      x = x & is.na(findOverlaps(query, excluderegions, select = "first"))
      frags = frags[x, ]
    }

    # Adjust frags$end to the Tn5 insertion site by subtracting Tn5 "width".
    frags$end = frags$end - 9

    # Count overlap at 5' end of each fragment.
    frags = frags %>%
      group_by(BC) %>%
      mutate(overlapcount = .overlapwithprecedingcount(start, end, TRUE)) %>%
      ungroup()

    # Summarize per BC.
    fragsbyBC = frags %>%
      group_by(BC) %>%
      summarize(nfrags = n(),
                depth1 = sum(overlapcount == 0),
                depth2 = sum(overlapcount == 1),
                depth3 = sum(overlapcount == 2),
                depth4 = sum(overlapcount == 3),
                depth5 = sum(overlapcount == 4),
                depth6 = sum(overlapcount == 5))

    sumoverlaplist = c(sumoverlaplist, list(fragsbyBC))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  sumoverlap =
    do.call(rbind, sumoverlaplist) %>%
    group_by(BC) %>%
    summarize(nfrags = sum(nfrags),
              depth1 = sum(depth1),
              depth2 = sum(depth2),
              depth3 = sum(depth3),
              depth4 = sum(depth4),
              depth5 = sum(depth5),
              depth6 = sum(depth6))
  sumoverlap = sumoverlap %>%
    rename(barcode = BC)
  return(sumoverlap)
}

# A utility function.
# The fragments must be sorted by start, end.
.overlapwithprecedingcount =
  function (start, end, include) {
    ct = rep(NA, length(start))
    if (length(include) == 1) {
      include = rep(include, length(start))
    }
    unfinishedends = c()
    for (i in 1:length(start)) {
      if (!is.na(include[i]) & include[i]) {
        unfinishedends = unfinishedends[unfinishedends >= start[i]]
        ct[i] = length(unfinishedends)
        unfinishedends = c(unfinishedends, end[i])
      }
    }
    return(ct)
  }
