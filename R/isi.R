# Methods to compute firing regularity based on Inter Spike Intervals (ISI)

##' Returns coefficient of local variation between consecutive ISIs, computed
##' with formula:
##' Lv = 3 / (n - 1) sum_{i=1}^{n-1} (ISI_i - ISI_{i+1}) ^ 2
##'     / (ISI_i + ISI_{i+1})^2
##'
##' Coefficient has value 0 when there is no variation between ISI, and is close
##' to 1 for irregular ISI.
##'
##' The measure proposed in:
##' Shinomoto S, Shima K, Tanji J (2003) Differences in spiking patterns among
##' cortical neurons. Neural Comput 15:2823–2842.
##' @param spike.train vector of spike times
##' @return coefficient of local variation of ISIs or 0 if not enough spikes
isi.local.variation <- function(spike.train) {
  isi <- diff(spike.train)
  n <- length(isi)
  if (n < 2) {
    return(0)
  }

  next.isi <- isi[2:n]
  isi <- isi[1:n-1]

  ret <- 3 / (n - 1) * sum((isi - next.isi)^2 / (isi + next.isi)^2)
  ret
}


##' Returns Spearman's rank correlation coefficient of ISI. The spike train is
##' split into chunks and the result is a mean of coefficients calculated per
##' each chunk.
##'
##' The coefficient measures monotonicity of ISI, and is calculated according
##' to formula:
##' rho = n / (n - 1) * sum_{i=1}^{i=n-1} (r_i - mean(r)) * (r_{i+1} - mean(r))
##'     / sum_{i=1}^{i=n} (r_i - mean(r)^2)
##' where r_i is a rank order of ISI i, equal ISI ranks are assigned average
##' rank
##'
##' @param spike.train vector of spike times
##' @param chunk.length length for chunks used for computation
##' @return Spearman's rank correlation coefficient or 0 if not enough spikes
isi.spearman.rank.corr <- function(spike.train, chunk.length=100) {
  isi <- diff(spike.train)
  n <- length(isi)

  chunk.count <- ceiling(n / chunk.length)
  if (n < 2) {
    return(0)
  }

  rho <- rep(0, chunk.count)
  for (i in 1:chunk.count) {
    index.start <- (i - 1) * chunk.length + 1
    index.end <- min(i * chunk.length, n)
    if (index.end - index.start < 2) {
      rho[i] <- 0
    } else {
      rho[i] <- cor(isi[index.start:(index.end - 1)],
                    isi[(index.start + 1):index.end],
                    method = "spearman")
    }
  }
  return(mean(rho))
}


##' Returns log(shape) and log(rate) of gamma distributions fitted to ISI of the
##' spike train. The spike train is split into chunks and gamma distribution for
##' each of them is estimated using maximum-likelihood. The result is a mean of
##' log of the estimated gamma distributions.
##'
##' The rate approximates mean firing rate, shape is a measure of firing
##' regularity:
##' - for firing in bursts: log(shape) < 1
##' - for Poisson distributed firing: log(shape) = 1
##' - for firing with peak of ISI: log(shape) > 1
##'
##' Calculation of the measure adopted from:
##' Y. Mochizuki et al., Similarity in neuronal firing regimes across mammalian
##' species. Journal of Neuroscience (2016) in press.
##'
##' @param spike.train vector of spike times
##' @param chunk.length length of chunks used for estimation
##' @return list with log(shape) and log(rate) with the means of fitted gamma
##'     distributions
isi.gamma <- function(spike.train, chunk.length=100) {
  require(MASS)
  isi <- diff(spike.train)
  n <- length(isi)

  chunk.count <- ceiling(n / chunk.length)
  if (chunk.count == 0) {
    return(list(logshape=NaN, lograte=NaN))
  }

  ret <- list(logshape=0, lograte=0)
  for (i in 1:chunk.count) {
    index.start <- (i - 1) * chunk.length + 1
    index.end <- min(i * chunk.length, n)
    if (index.end - index.start >= 2) {
      estimate <- fitdistr(isi[index.start:index.end], "gamma", lower=0.001)$estimate
      ret$logshape = ret$logshape + log(estimate[1]) / chunk.count
      ret$lograte = ret$lograte + log(estimate[2]) / chunk.count
    }
  }

  ret
}

