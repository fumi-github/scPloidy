# Core functions to infer ploidy

if (getRversion() >= "2.15.1") {
  # For nimble
  utils::globalVariables(c("alpha1", "averagedepth1", "equals", "ind", "ploidylevels"))
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

### EM ALGORITHM FOR MIXTURES
# We superficially (and possibly robustly) model
# as mixtures of multinomial distributions
# We first try with (depth2, depth3, depth4),
# but if the clusters don't separate,
# next try (depth3, depth4, depth5), and so on.
inferpem = function (fragmentoverlap, levels, s, epsilon, subsamplesize) {
  if (ncol(fragmentoverlap) < 6) { return(NA) }
  for (j in 4:(ncol(fragmentoverlap) - 2)) {
    fragmentoverlapsubmatrix =
      as.matrix(fragmentoverlap[, 0:2 + j])
    lambda = NULL
    theta = NULL
    if (is.numeric(subsamplesize)) {
      while (subsamplesize < nrow(fragmentoverlapsubmatrix)) {
        set.seed(s)
        em.out.small = multmixEM(
          y = fragmentoverlapsubmatrix[
            sample(nrow(fragmentoverlapsubmatrix), subsamplesize), ],
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
      y = fragmentoverlapsubmatrix,
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
  return(p.em)
}

### K-MEANS POST-PROCESSING OF MOMENT
inferpkmeans = function (fragmentoverlap, levels, p.moment) {
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
}

### BAYESIAN
# TODO cells with no observation (rowSums(data) == 0) might cause error.
inferpbayes = function (data, levels, prop, inits) {

  Code = nimbleCode({
    alpha1 ~ dbeta(shape1 = 0.5, shape2 = 0.5)
    for (i in 1:Ncell) {
      ind[i] ~ dcat(ploidyprior[])
      ploidy[i] <- ploidylevels[ind[i]]
      for (j in 1:6) {
        probbinom1raw[i, j] <- dbinom(x = j, prob = alpha1, size = ploidy[i], log = 0)
      }

      # avoid if
      probpois1raw[i, 1]  <- dpois(x = 1, lambda = averagedepth1[i], log = 0) + equals(averagedepth1[i], 0)
      for (j in 2:6) {
        probpois1raw[i, j]  <- dpois(x = j, lambda = averagedepth1[i], log = 0)
      }

      probbinom1sum[i] <- sum(probbinom1raw[i, 1:6])
      probpois1sum[i]  <- sum(probpois1raw[i, 1:6])
      for (j in 1:6) {
        prob1[i, j] <-
          prop * (probbinom1raw[i, j] / probbinom1sum[i]) +
          (1 - prop) * (probpois1raw[i, j] / probpois1sum[i])
      }
      y[i, 1:6] ~ dmulti(prob = prob1[i, 1:6], size = Nfrag1[i])
    }
  })

  Ncell = nrow(data)
  T1_1 = as.numeric(data[, 1:6] %*% seq(1, 6))
  T2_1 = as.numeric(data[, 1:6] %*% (seq(1, 6)^2))
  T2T1_1 = T2_1 / T1_1 - 1
  Consts = list(
    prop = prop,
    ploidyprior = rep(1/length(levels), length(levels)),
    Ncell = Ncell,
    Nfrag1 = rowSums(data[, 1:6]),
    averagedepth1 = T2T1_1)
  Data = list(
    ploidylevels = levels,
    y = data)
  Inits = list(
    ind = match(inits, levels))
  mcmc.out = nimbleMCMC(
    code = Code,
    constants = Consts,
    data = Data,
    inits = Inits,
    nchains = 4, niter = 30000, nburnin = 20000,
    setSeed = TRUE,
    summary = TRUE, WAIC = TRUE,
    monitors = c('alpha1', 'ind'))

  # v = "alpha1"
  # ggplot(data = data.frame(
  #   x = 1:nrow(mcmc.out$samples$chain1),
  #   y1 = mcmc.out$samples$chain1[, v],
  #   y2 = mcmc.out$samples$chain2[, v],
  #   y3 = mcmc.out$samples$chain3[, v],
  #   y4 = mcmc.out$samples$chain4[, v]),
  #   aes(x = x)) +
  #   geom_line(aes(y = y1)) +
  #   geom_line(aes(y = y2)) +
  #   geom_line(aes(y = y3)) +
  #   geom_line(aes(y = y4))
  # mcmc.out$summary$all.chains
  # plogis(mcmc.out$summary$all.chains["alpha1", "Mean"])

  alpha1 = (mcmc.out$summary$all.chains)["alpha1", "Median"]
  x = grep("^ind", rownames(mcmc.out$summary$all.chains))
  ploidy.bayes = levels[round((mcmc.out$summary$all.chains)[x, "Median"])]
  return(list(
    ploidy.bayes = ploidy.bayes,
    alpha1 = alpha1))

}
