#' Infer Copy Number Variations (CNVs) in Cancer Cells from ATAC-seq Fragment Overlap
#'
#' @param fragmentoverlap Frequency of fragment overlap in each cell-window
#' computed by the function \code{fragmentoverlapcount}.
#' \code{barcode} should be named as \code{AAACGAAAGATTGACA-1.window_1},
#' which represents cell \code{AAACGAAAGATTGACA-1} and window \code{window_1}.
#' The format is "cell barcode", ".window_" and integer.
#' @param windowcovariates Chromosomal windows for which copy number
#' gain/loss are initially inferred. Required columns are chr, start, end,
#' window (for example, \code{window_1}) and peaks.
#' Peaks is a numeric column representing chromatin accessibility.
#' @param levels Possible values of ploidy. For example,
#' \code{c(2, 4)} if the cells can be diploids or tetraploids.
#' The values must be larger than one.
#' @param nfragspercellmin Minimum number of fragments for a cell-window to be eligible.
#' @param nfragspercellmax Maximum number of fragments for a cell-window to be eligible.
#' @param deltaBICthreshold Only the CNVs with deltaBIC smaller than this threshold
#' are adopted.
#' @return A list with two elements.
#' \code{CNV} is a data frame of the CNVs identified in the dataset.
#' \code{cellwindowCN} is a data frame indicating the ploidy for each cell
#' and the inferred standardized copy number for each cell-window.
#'
#' @importFrom dplyr filter group_by mutate n summarize ungroup
#' @importFrom MASS rlm
#' @importFrom rlang warn
#' @importFrom stats median
#' @importFrom tibble tibble
#' @export
cnv = function(fragmentoverlap,
               windowcovariates,
               levels = c(2, 4),
               nfragspercellmin = 5000,
               nfragspercellmax = 10^5.5,
               deltaBICthreshold = 0) {
  spanlist = 2:10

  fragmentoverlap = fragmentoverlap[, 1:8]
  fragmentoverlap$cell =
    sub(".window.*", "", fragmentoverlap$barcode)
  fragmentoverlap$window =
    sub(".*window", "window", fragmentoverlap$barcode)
  fragmentoverlap$w =
    as.numeric(sub("window_", "", fragmentoverlap$window))
  windowcovariates$w =
    as.numeric(sub("window_", "", windowcovariates$window))

  # QC of cells
  x = fragmentoverlap %>%
    group_by(.data$cell) %>%
    summarize(n = n(), s = sum(.data$nfrags))
  fragmentoverlap =
    fragmentoverlap %>%
    group_by(.data$cell) %>%
    mutate(n = n(), s = sum(.data$nfrags)) %>%
    filter(.data$n > max(x$n) * 0.85) %>%
    filter(.data$s >= nfragspercellmin) %>%
    filter(.data$s <= nfragspercellmax) %>%
    mutate(n = NULL, s = NULL) %>%
    ungroup()

  # Infer ploidy of cells
  fragmentoverlapcell =
    fragmentoverlap %>%
    group_by(.data$cell) %>%
    summarize(nfrags = sum(.data$nfrags),
              depth1 = sum(.data$depth1),
              depth2 = sum(.data$depth2),
              depth3 = sum(.data$depth3),
              depth4 = sum(.data$depth4),
              depth5 = sum(.data$depth5),
              depth6 = sum(.data$depth6))
  x = as.matrix(fragmentoverlapcell[, 3:8])
  T1 = as.numeric(x %*% seq(1, ncol(x)))
  T2 = as.numeric(x %*% (seq(1, ncol(x))^2))
  logT2T1 = log(T2 / T1 - 1)
  logT2T1capped = .cap(logT2T1)
  x = inferpmoment(logT2T1capped, levels)
  fragmentoverlap$p.moment.cell =
    x$p.moment[match(fragmentoverlap$cell, fragmentoverlapcell$cell)]
  rm(T1, T2, logT2T1, logT2T1capped, x)

  # QC of windows
  x = fragmentoverlap %>%
    group_by(.data$window) %>%
    summarize(n = n())
  fragmentoverlap =
    fragmentoverlap %>%
    group_by(.data$window) %>%
    mutate(n = n()) %>%
    filter(.data$n > max(x$n) * 0.85) %>%
    mutate(n = NULL) %>%
    ungroup()

  # Compute T2T1 per cell-window
  x = as.matrix(fragmentoverlap[, 3:8])
  T1 = as.numeric(x %*% seq(1, ncol(x)))
  T2 = as.numeric(x %*% (seq(1, ncol(x))^2))
  fragmentoverlap$T2T1 = T2 / T1 - 1

  # Further QC of cells
  fragmentoverlap =
    fragmentoverlap %>%
    group_by(.data$cell) %>%
    mutate(T2T1median = median(.data$T2T1)) %>%
    filter(.data$T2T1median > 0) %>%
    mutate(T2T1median = NULL) %>%
    ungroup()

  # Attach window features
  x = match(fragmentoverlap$window, windowcovariates$window)
  fragmentoverlap$chr   = windowcovariates$chr[x]
  fragmentoverlap$peaks = windowcovariates$peaks[x]

  # Correct T2T1 by ATAC-seq accessibility per window
  # Correction is based on lowest ploidy cells.
  # p: ploidy, s: probability of sequencing
  # T2T1              = (p_cell_window - 1) * s_window
  # fitted.values     = (p_minlevel    - 1) * s_window
  # T2T1scaled        = { (p_cell_window - 1) * s_window } / { (p_minlevel - 1) * s_window }
  #                   = { (p_cell_window - 1)         } / { (p_minlevel - 1)  }
  # T2T1scaledmedian  = { (p_cell        - 1)         } / { (p_minlevel - 1)  }
  # pwindownormalized = p_minlevel       * (p_cell_window ) / (p_cell )     ; intuitive
  # pwindownormalized = (p_minlevel - 1) * (p_cell_window ) / (p_cell ) + 1 ; stabilized
  x = fragmentoverlap[
    fragmentoverlap$p.moment.cell == min(levels), ]
  a0 = rlm(T2T1 ~ peaks, data = x)
  if (a0$coefficients["peaks"] <= 0 | min(a0$fitted.values) <= 0) {
    warn("Inappropriate fitting of T2T1 ~ peaks + intercept. Intercept set to zero.")
    a0 = rlm(T2T1 ~ peaks +  0, data = x)
  }
  x$fitted.values = a0$fitted.values
  x = x %>%
    dplyr::select(.data$window, .data$fitted.values) %>%
    dplyr::distinct()
  fragmentoverlap$fitted.values =
    x$fitted.values[match(fragmentoverlap$window, x$window)]
  fragmentoverlap = fragmentoverlap[!is.na(fragmentoverlap$fitted.values), ]
  # TODO improve for cases min(levels) != 2
  fragmentoverlap$T2T1scaled = fragmentoverlap$T2T1 / fragmentoverlap$fitted.values
  # Normalize cell-window T2T1 by the value per cell and convert to ploidy
  fragmentoverlap =
    fragmentoverlap %>%
    group_by(.data$cell) %>%
    mutate(T2T1scaledmedian = median(.data$T2T1scaled)) %>%
    # +1 in numerator and denominator is practically for stabilization. Theoretically inaccurate.
    mutate(pwindownormalized = min(levels) * (.data$T2T1scaled + 1) / (.data$T2T1scaledmedian + 1)) %>%
    ungroup()

  ### Detect CNV
  # index and vec are vectors of same length
  get_values_in_span = function(index, vec, span) {
    lapply(index, function(x) { vec[index >= x & index <= x + span - 1] })
  }

  # A. Fit different probabilities of variant occurrence in
  #   (windows in span) x (across up to n-th cell) vs all other cell-windows
  # B. Fit single probability of variant occurrence in all cell-windows
  # Minimize: deltaBIC = BIC_A - BIC_B
  mindeltaBIC = function(celllist, variantlist, totalvariant, totalnotvariant) {
    # cl = celllist[[1]][[1]]
    # vl = variantlist[[1]][[1]]
    cl = celllist[[1]]
    vl = variantlist[[1]]
    o = order(unlist(lapply(vl, mean)), decreasing = TRUE)
    cl = cl[o]
    vl = vl[o]
    a = unlist(lapply(vl, function(x) { sum(x == 1) }))
    a = cumsum(a) # up to n-th cell
    b = unlist(lapply(vl, function(x) { sum(x == 0) }))
    b = cumsum(b) # up to n-th cell
    c = totalvariant - a
    d = totalnotvariant - b
    deltaloglikelihood =
      ifelse(a > 0, a * log(a / (a + b)), 0) +
      ifelse(b > 0, b * log(b / (a + b)), 0) +
      ifelse(c > 0, c * log(c / (c + d)), 0) +
      ifelse(d > 0, d * log(d / (c + d)), 0) -
      totalvariant * log(totalvariant / (totalvariant + totalnotvariant)) -
      totalnotvariant * log(totalnotvariant / (totalvariant + totalnotvariant))
    deltaBIC =
      log(totalvariant + totalnotvariant) - 2 * deltaloglikelihood
    i = which.min(deltaBIC)
    return(list(deltaBIC[i], unlist(cl[1:i])))
  }

  cnv = function(fragmentoverlap, spanlist) {
    totalvariant = sum(fragmentoverlap$variant == 1)
    totalnotvariant = sum(fragmentoverlap$variant == 0)
    result = tibble()
    for (span in spanlist) {
      print(paste0("Computing span = ", span))
      wvariant = fragmentoverlap %>%
        dplyr::select(.data$cell, .data$w, .data$chr, .data$variant) %>%
        group_by(.data$cell, .data$chr) %>%
        filter(.data$w <= max(.data$w) - span + 1) %>% # omit w at end of chromosome
        mutate(variant_in_span = get_values_in_span(.data$w, .data$variant, span)) %>%
        ungroup() %>%
        group_by(.data$w) %>%
        summarize(celllist = list(.data$cell), variantlist = list(.data$variant_in_span))
      wvariant$deltaBIC = NA
      wvariant$deltaBICcelllist = NA
      for (i in 1:nrow(wvariant)) {
        x = mindeltaBIC(wvariant$celllist[i], wvariant$variantlist[i], totalvariant, totalnotvariant)
        wvariant$deltaBIC[i] = x[[1]]
        wvariant$deltaBICcelllist[i] = list(x[[2]])
      }
      wvariant$span = span
      result = rbind(result, wvariant)
    }
    return(result)
  }

  # copy number gain variant
  fragmentoverlap$variant = 1 * (fragmentoverlap$pwindownormalized >= min(levels) + 1)
  resultgain = cnv(fragmentoverlap, spanlist)
  # copy number loss variant
  fragmentoverlap$variant = 1 * (fragmentoverlap$pwindownormalized <= min(levels) - 1)
  resultloss = cnv(fragmentoverlap, spanlist)

  x = resultgain; x$type = "gain"
  y = resultloss; y$type = "loss"
  result =
    rbind(x, y) %>%
    dplyr::arrange(.data$deltaBIC)
  resultnonoverlapping = tibble()
  while (nrow(result) > 0) {
    resultnonoverlapping = rbind(
      resultnonoverlapping,
      result[1, ])
    w = result$w[1]
    span = result$span[1]
    result = result[result$w + result$span - 1 < w | result$w > w + span - 1, ]
  }

  resulttop =
    resultnonoverlapping %>%
    dplyr::arrange(.data$w) %>%
    dplyr::select(.data$w, .data$deltaBIC, .data$deltaBICcelllist, .data$span, .data$type) %>%
    filter(.data$deltaBIC < deltaBICthreshold)
  x = match(resulttop$w, windowcovariates$w)
  y = match(resulttop$w + resulttop$span - 1, windowcovariates$w)
  resulttop$chr = windowcovariates$chr[x]
  resulttop$start = windowcovariates$start[x]
  resulttop$end   = windowcovariates$end[y]

  fragmentoverlap$cnv = 0
  for (i in 1:nrow(resulttop)) {
    x =
      fragmentoverlap$w >= resulttop$w[i] &
      fragmentoverlap$w <= resulttop$w[i] + resulttop$span[i] - 1
    fragmentoverlap$cnv[x] =
      ifelse(resulttop$type[i] == "gain", 1, -1)
  }
  fragmentoverlap$pwindownormalizedcleaned = min(levels)
  x = (fragmentoverlap$cnv == 1)
  fragmentoverlap$pwindownormalizedcleaned[x] =
    pmax(fragmentoverlap$pwindownormalized[x], min(levels))
  x = (fragmentoverlap$cnv == -1)
  fragmentoverlap$pwindownormalizedcleaned[x] =
    pmin(fragmentoverlap$pwindownormalized[x], min(levels))

  resulttop = as.data.frame(resulttop[, c("chr", "start", "end", "type", "deltaBIC")])
  x = as.data.frame(fragmentoverlap[, c("barcode", "p.moment.cell", "pwindownormalizedcleaned")])
  colnames(x)[2:3] = c("ploidy.moment.cell", "CN")
  return(list(CNV = resulttop,
              cellwindowCN = x))
}
