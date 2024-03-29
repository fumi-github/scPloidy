#' Count Overlap of ATAC-seq Fragments
#'
#' @param file Filename of the file for ATAC-seq fragments.
#' The file must be block gzipped (using the \code{bgzip} command)
#' and accompanied with the index file (made using the \code{tabix} command).
#' The uncompressed file must be a tab delimited file,
#' where each row represents one fragment.
#' The first four columns are chromosome name, start position, end position,
#' and barcode (i.e., name) of the cell including the fragment.
#' The remaining columns are ignored.
#' See vignette for details.
#' @param targetregions GRanges object for the regions where overlaps are counted.
#' Usually all of the autosomes.
#' If there is memory problem, split a chromosome into smaller chunks,
#' for example by 10 Mb.
#' The function loads each element of \code{targetregions} sequentially,
#' and smaller elements require less memory.
#' @param excluderegions GRanges object for the regions to be excluded.
#' Simple repeats in the genome should be listed here,
#' because repeats can cause false overlaps.
#' A fragment is discarded if its 5' or 3' end is located in \code{excluderegions}.
#' If \code{NULL}, fragments are not excluded by this criterion.
#' @param targetbarcodes Character vector for the barcodes of cells to be analyzed,
#' such as those passing quality control.
#' If \code{NULL}, all barcodes in the input file are analyzed.
#' @param Tn5offset Numeric vector of length two.
#' The enzyme for ATAC-seq is a homodimer of Tn5.
#' The transposition sites of two Tn5 proteins are 9 bp apart,
#' and the (representative) site of accessibility is in between.
#' If the start and end position of your input file is taken from BAM file,
#' set the paramater to \code{c(4, -5)} to adjust the offset.
#' Alternatively, values such as \code{c(0, -9)} could generate similar results;
#' what matters the most is the difference between the two numbers.
#' The fragments.tsv.gz file generated by 10x Cell Ranger already adjusts the shift
#' but is recorded as a BED file. In this case, use \code{c(1, 0)} (default value).
#' If unsure, set to \code{"guess"},
#' in which case the program returns a guess.
#' @return A tibble with each row corresponding to a cell.
#' For each cell, its barcode, the total count of the fragments \code{nfrag},
#' and the count distinguished by overlap depth are given.
#'
#' @importFrom dplyr arrange desc distinct filter group_by lag left_join mutate n rename summarize ungroup
#' @importFrom GenomicRanges findOverlaps makeGRangesFromDataFrame
#' @importFrom rlang .data
#' @importFrom tidyr pivot_wider
#' @importFrom Rsamtools scanTabix TabixFile
#' @importFrom stats ksmooth
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar
#' @export
fragmentoverlapcount = function (file,
                                 targetregions,
                                 excluderegions = NULL,
                                 targetbarcodes = NULL,
                                 Tn5offset = c(1, 0)) {
  sumoverlaplist = list()
  tbx = TabixFile(file = file)
  pb = txtProgressBar(max = length(targetregions), style = 3)
  for (i in 1:length(targetregions)) {

    # Load fragments.
    res = scanTabix(tbx, param = targetregions[i])
    frags = read.csv(textConnection(res[[1]]),
                     sep = "\t",
                     header = FALSE)
    frags = frags[, 1:4]
    colnames(frags) = c("chr", "start", "end", "BC")
    rm(res)

    if (! is.null(targetbarcodes)) {
      frags = frags[frags$BC %in% targetbarcodes, ]
    }
    if (nrow(frags) == 0) { next() }

    # Discard "semi-duplicate" fragments;
    # Hypothesized Tn5 transposition only at one strand.
    frags = frags %>%
      group_by(.data$BC) %>%
      # Between chrX:100-200 and chrX:100-300, only retain first
      arrange(.data$start, .data$end) %>%
      distinct(.data$start, .keep_all = TRUE) %>%
      # Between chrX:100-300 and chrX:200-300, only retain last
      arrange(desc(.data$end), desc(.data$start)) %>%
      distinct(.data$end, .keep_all = TRUE) %>%
      ungroup() %>%
      arrange(.data$start, .data$end, .data$BC)

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
      if (nrow(frags) == 0) { next() }
    }

    # Adjust Tn5 site offset
    if (identical(Tn5offset, "guess")) {
      x = frags %>%
        group_by(.data$BC) %>%
        mutate(overlap = (lag(.data$end) - .data$start + 1)) %>%
        ungroup()
      x = x$overlap
      x = x[abs(x) <= 18]
      if (length(x) < 20) {
        stop('Error: datasize is too small for guessing Tn5offset')
      }
      x = table(x)
      x = as.numeric(names(x)[which.max(x)])
      x = c(0, -x)
      x = x - round(mean(x))
      setTxtProgressBar(pb, length(targetregions))
      close(pb)
      return(x)
    } else {
      frags$start = frags$start + Tn5offset[1]
      frags$end   = frags$end   + Tn5offset[2]
    }

    # Count overlap at 5' end of each fragment.
    frags = frags %>%
      group_by(.data$BC) %>%
      mutate(overlapcount = .overlapwithprecedingcount(.data$start, .data$end, TRUE)) %>%
      ungroup()

    # Summarize per BC.
    fragsbyBC = frags %>%
      group_by(.data$BC) %>%
      summarize(nfrags = n(),
                depth1 = sum(.data$overlapcount == 0),
                depth2 = sum(.data$overlapcount == 1),
                depth3 = sum(.data$overlapcount == 2),
                depth4 = sum(.data$overlapcount == 3),
                depth5 = sum(.data$overlapcount == 4),
                depth6 = sum(.data$overlapcount == 5))

    # Compute smoothed distance to the next fragment (irrelevant to BC),
    # which is the inverse of chromatin accessibility.
    compute_smoothed_distance <-
      function(frags) {
        smoothed_starts <- ksmooth(1:nrow(frags), frags$start, bandwidth = 200)$y
        differences <- diff(smoothed_starts)
        differences <- c(differences, differences[length(differences)])
        frags$bptonext <- pmax(differences, 0)
        return(frags)
      }

    # Convert bptonext into a list format
    bptonext_as_list <-
      function(frags) {
        list_data <- frags %>%
          mutate(depth = .data$overlapcount + 1) %>%
          dplyr::filter(.data$depth <= 6) %>%
          group_by(.data$BC, .data$depth) %>%
          summarize(bptonext = list(.data$bptonext), .groups = "drop")
        widened_data <- pivot_wider(
          list_data,
          names_prefix = "bptonextdepth",
          names_from = .data$depth,
          values_from = .data$bptonext
        )
        for (i in setdiff(paste0("bptonextdepth", 1:6), colnames(widened_data))) {
          widened_data[[i]] <- list(NULL)
        }
        return(widened_data)
      }

    frags <- compute_smoothed_distance(frags)
    bptonext_list_data <- bptonext_as_list(frags)
    fragsbyBC <-
      left_join(fragsbyBC, bptonext_list_data, by = 'BC')

    sumoverlaplist = c(sumoverlaplist, list(fragsbyBC))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  if (identical(sumoverlaplist, list())) {
    stop('Error: no fragments remained after filtering')
  }

  sumoverlap =
    do.call(rbind, sumoverlaplist) %>%
    group_by(.data$BC) %>%
    summarize(nfrags = sum(.data$nfrags),
              depth1 = sum(.data$depth1),
              depth2 = sum(.data$depth2),
              depth3 = sum(.data$depth3),
              depth4 = sum(.data$depth4),
              depth5 = sum(.data$depth5),
              depth6 = sum(.data$depth6),
              bptonextdepth1 = list(do.call(c, .data$bptonextdepth1)),
              bptonextdepth2 = list(do.call(c, .data$bptonextdepth2)),
              bptonextdepth3 = list(do.call(c, .data$bptonextdepth3)),
              bptonextdepth4 = list(do.call(c, .data$bptonextdepth4)),
              bptonextdepth5 = list(do.call(c, .data$bptonextdepth5)),
              bptonextdepth6 = list(do.call(c, .data$bptonextdepth6)))
  sumoverlap = sumoverlap %>%
    rename(barcode = .data$BC)
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

#' Infer Ploidy from ATAC-seq Fragment Overlap
#'
#' @param fragmentoverlap Frequency of fragment overlap in each cell
#' computed by the function \code{fragmentoverlapcount}.
#' @param levels Possible values of ploidy. For example,
#' \code{c(2, 4)} if the cells can be diploids or tetraploids.
#' The values must be larger than one.
#' @param s Seed for random numbers used in EM algorithm.
#' @param epsilon Convergence criterion for the EM algorithm.
#' @param subsamplesize EM algorithm becomes difficult to converge
#' when the number of cells is very large.
#' By setting the parameter (e.g. to 1e4),
#' we can run EM algorithm iteratively,
#' first for \code{subsamplesize} randomly sampled cells,
#' next for twice the number of cells in repetition.
#' The inferred lambda/theta parameters are used as the initial value
#' in the next repetition.
#' @return A data.frame with each row corresponding to a cell.
#' For each cell, its barcode, ploidy inferred by moment method,
#' the same with additional K-means clustering,
#' and ploidy inferred by EM algorithm of mixture are given.
#' I recommend using \code{ploidy.moment}.
#'
#' @importFrom matrixStats colMins
#' @importFrom mixtools multmixEM
#' @importFrom stats kmeans optimize quantile
#' @export
ploidy = function (fragmentoverlap,
                   levels,
                   s = 100,
                   epsilon = 1e-08,
                   subsamplesize = NULL) {
  if (min(levels) <= 1) {
    stop('Error: elements of levels must be larger than one')
  }

  ### MOMENT BASED METHOD
  # We model the overlapping of fragments by binomial distribution:
  # ------------------------------------------------------------
  # binomial distribution   | overlap of fragments   | parameter
  # ------------------------------------------------------------
  # one observation         | 5' end of a fragment   |
  # number of trials (size) | ploidy                 | p
  # number of success       | depth of overlap       |
  # probability of success  | probability of overlap | s
  # ------------------------------------------------------------
  # Under a predetermined p, for each cell, we estimate s based on
  # the overlap depth observed in the fragments belonging to the cell.
  # Since we cannot properly count observations with zero success,
  # we model as truncated binomial distribution.
  # We use the moment method in Paul R. Rider (1955).

  # Computes a "capped" version of the log-transformed T2T1 values.
  # The capping is done based on the interquartile range of the log-transformed values.
  cap = function (logT2T1) {
    x = 2 * quantile(logT2T1, 0.75) -
      quantile(logT2T1, 0.5)
    logT2T1capped = pmin(logT2T1, x)
    x = 2 * quantile(logT2T1, 0.25) -
      quantile(logT2T1, 0.5)
    logT2T1capped = pmax(logT2T1capped, x)
    # Although unlikely, adjust if -Inf remains.
    x = min(logT2T1capped[is.finite(logT2T1capped)])
    logT2T1capped = pmax(logT2T1capped, x)
    return(logT2T1capped)
  }

  inferpmoment = function (logT2T1capped, levels) {
    m = matrix(
      logT2T1capped,
      nrow = length(levels),
      ncol = length(logT2T1capped),
      byrow = TRUE)
    # Optimizes an offset for the log-transformed T2T1 values so that the
    # squared differences from the log(levels-1) are minimized.
    offsetoptimize =
      optimize(
        function (o) {
          x = m + o - log(levels - 1)
          return(sum(matrixStats::colMins(x^2))) },
        lower = min(log(levels - 1)) -
          max(logT2T1capped),
        upper = max(log(levels - 1)) -
          min(logT2T1capped))
    # Infers the ploidy using the optimized offset and the provided levels.
    # It computes the closest level for each log-transformed value by
    # minimizing the absolute differences.
    p.moment =
      apply(
        abs(m + offsetoptimize$minimum - log(levels - 1)),
        2,
        which.min)
    p.moment = levels[p.moment]
    return(list(
      p.moment = p.moment,
      offset = offsetoptimize$minimum))
  }

  x = as.matrix(fragmentoverlap[, 3:8])
  T1 = as.numeric(x %*% seq(1, ncol(x)))
  T2 = as.numeric(x %*% (seq(1, ncol(x))^2))
  logT2T1 = log(T2 / T1 - 1)
  logT2T1capped = cap(logT2T1)
  x = inferpmoment(logT2T1capped, levels)
  p.moment = x$p.moment
  offset = x$offset

  # exp(offset) is the estimate for 1/s
  p.momentfractional =
    exp(logT2T1) * exp(offset) + 1

  ### EM ALGORITHM FOR MIXTURES
  # We superficially (and possibly robustly) model
  # as mixtures of multinomial distributions
  # We first try with (depth2, depth3, depth4),
  # but if the clusters don't separate,
  # next try (depth3, depth4, depth5), and so on.
  for (j in 4:6) {
    sumoverlapsubmatrix =
      as.matrix(fragmentoverlap[, 0:2 + j])
    lambda = NULL
    theta = NULL
    if (is.numeric(subsamplesize)) {
      while (subsamplesize < nrow(sumoverlapsubmatrix)) {
        set.seed(s)
        em.out.small = multmixEM(
          y = sumoverlapsubmatrix[
            sample(nrow(sumoverlapsubmatrix), subsamplesize), ],
          lambda = lambda,
          theta = theta,
          k = length(levels),
          epsilon = epsilon)
        lambda = em.out.small$lambda
        theta = em.out.small$theta
        subsamplesize = 2 * subsamplesize
      }
    }
    set.seed(s)
    em.out = multmixEM(
      y = sumoverlapsubmatrix,
      lambda = lambda,
      theta = theta,
      k = length(levels),
      epsilon = epsilon)
    if (max(em.out$lambda) < 0.99) { break }
  }
  p.em = apply(em.out$posterior, 1, which.max)
  # EM is simple clustering and unaware of the labeling in levels.
  # We infer the labeling from the last element of theta,
  # which represents the overlaps of largest depth used for clustering.
  p.em = ( sort(levels)[ rank(em.out$theta[, 3]) ] )[p.em]

  ### K-MEANS POST-PROCESSING OF MOMENT
  x = log10(
    (fragmentoverlap[, 4:6] + 1) / fragmentoverlap$nfrags)
  kmclust =
    kmeans(
      x,
      do.call(
        rbind,
        tapply(
          as.list(as.data.frame(t(x))),
          p.moment,
          function (x) {rowMeans(do.call(cbind, x))})))
  p.kmeans = levels[kmclust$cluster]

  return(data.frame(
    barcode = fragmentoverlap$barcode,
    ploidy.moment = p.moment,
    ploidy.momentfractional = p.momentfractional,
    ploidy.kmeans = p.kmeans,
    ploidy.em = p.em))
}
