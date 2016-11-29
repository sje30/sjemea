################################################################################
# Tests for coefficient of ISI Loval Variation                                 #
################################################################################
context("Test coefficient of ISI Loval Variation when no spikes")
expect_equal(isi.local.variation(c()), 0)

context("Test coefficient of ISI Loval Variation when single spike")
expect_equal(isi.local.variation(c(1)), 0)

context("Test coefficient of ISI Loval Variation when single ISI")
expect_equal(isi.local.variation(c(1, 3)), 0)

context("Test coefficient of ISI Loval Variation when ISI constant")
expect_equal(isi.local.variation(1:10), 0)

context("Test coefficient of ISI Loval Variation when ISI vary")
expect_equal(isi.local.variation(c(0, 2, 3)), 1/3)
expect_equal(isi.local.variation(c(0, 2, 3, 5)), 1/3)
expect_equal(isi.local.variation(c(0, 2, 3, 5, 12)), 43/81)

################################################################################
# Tests for Spearman's rank correlation coefficient                            #
################################################################################
context("Test Spearman's rank correlation coefficient when no spikes")
expect_equal(isi.spearman.rank.corr(c()), 0)

context("Test Spearman's rank correlation coefficient when single spike")
expect_equal(isi.spearman.rank.corr(c(1)), 0)

context("Test Spearman's rank correlation coefficient when single ISI")
expect_equal(isi.spearman.rank.corr(c(1, 2)), 0)

context("Test Spearman's rank correlation coefficient when ISI monothically
        increasing")
expect_equal(isi.spearman.rank.corr(c(1, 2, 4, 7)), 1)

context("Test Spearman's rank correlation coefficient when ISI monothically
        decreasing")
expect_equal(isi.spearman.rank.corr(c(7, 4, 2, 1)), 1)

context("Test Spearman's rank correlation coefficient when ISI vary")
expect_equal(isi.spearman.rank.corr(c(0, 3, 4, 8, 10)), -0.5)

context("Test Spearman's rank correlation coefficient when averaged from two
        chunks")
rank.corr <- isi.spearman.rank.corr(c(0, 3, 4, 8, 10, 11, 13, 16),
                                    chunk.length = 4)
expect_equal(rank.corr, 0.25)

################################################################################
# Tests for gamma distribution for ISI                                         #
################################################################################
assert.gamma.estimate <- function(actual.estimate, expected.estimate) {
  expect_true(abs(actual.estimate$logshape - expected.estimate$logshape) < 0.2)
  expect_true(abs(actual.estimate$lograte - expected.estimate$lograte) < 0.1)
}

context("Test gamma distribution estimate for empty spike train")
e <- isi.gamma(c())
expect_equal(e$logshape, NaN)
expect_equal(e$lograte, NaN)

context("Test gamma distribution estimate for signle spike")
e <- isi.gamma(c(1))
expect_equal(e$logshape, NaN)
expect_equal(e$lograte, NaN)

context("Test gamma distribution fitted to ISI of entire spike train")
x <- rgamma(2000, shape = 5, rate = 0.1)
e <- isi.gamma(cumsum(x), chunk.length = 2000)
assert.gamma.estimate(e, list(logshape=log(5), lograte=log(0.1)))

context("Test gamma distribution fitted to ISI of entire spike train")
x <- c(rgamma(2000, shape = 5, rate = 0.1), rgamma(2000, shape=3, rate = 2))
e <- isi.gamma(cumsum(x), chunk.length = 2000)
expected.estimate = list(logshape=(log(5) + log(3)) / 2, 
                         lograte=(log(0.1) + log(2)) / 2)
assert.gamma.estimate(e, expected.estimate)
