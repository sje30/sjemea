## Test the tiling coefficients.

## if one train is empty, we get NaN.
expect_equal( tiling.corr(1:5, double(0)), NaN)

context("Tiling coefficent should return 1 for autocorrelated trains.")

poisson.train <- function(n=1000, rate=1, beg=0) {
  ## Generate a Poisson spike train with N spikes and firing rate RATE.
  ## BEG is time of start of recording
  ## Check that the histogram looks exponentially distributed.
  ## hist( diff(poisson.train()))
  x <- runif(n)
  isi <- log(x) / rate
  spikes <- beg - cumsum(isi)
  spikes
}

## Auto-correlations should be 1.

## Poisson spike train generated using the rule: t_i+1 = t_i -
## ln(x)/r, where x is drawn from Uniform distribuition and r is the
## rate.

## So in each case lets draw n=200 spikes from varying a rates.
n <-  5000
rates <- c(0.3, 1, 2, 5, 10)
for (r in rates) {
  t <- poisson.train(n, r, beg=300)
  expect_equal( tiling.corr(t, t), 1)
}

context("Tiling coefficent should return 0 for two independent Poisson trains.")

## We look at the distribution of many trials.  the tails of the
## distribution should still be close to zero.  Here we check that the
## 5% and 95% bins are less than 0.02 away from zero, and that they
## are opposite signs, so mirrored around zero.

tiling.ind <- function()  {
  n <- 3000; r <- 0.2
  ## Compute tiling for a pair of Poisson trains -- should be close to zero.
  a <- poisson.train(n, r, beg=300)
  b <- poisson.train(n, r, beg=300)
  tiling.corr(a, b)
}

coefs <- replicate(1000, tiling.ind())
hist(coefs)
percentiles <- quantile(coefs, probs=c(0.05, 0.95))
expect_true(max(abs(percentiles)) < 0.02)
expect_true( prod(percentiles) < 0)     #check opposite signs.


context("Introducing some correlation between pairs of trains.")

tiling.shared <- function(p.shared=0.6, n=2000, r=1) {
  master <- poisson.train(n, r)
  p.own <- (1-p.shared)/2
  p <- c( p.own, p.own, p.shared)
  ## Each spike is in one of three states with prob given by P vector above.
  ## 1: in train 1 only.
  ## 2: in train 2 only.
  ## 3: in both trains
  state <- sample(1:3, n, replace=TRUE, prob=p)
  a <- master[state != 2]
  b <- master[state != 1]
  tiling.corr(a, b)
}

p.shared <- seq(from=0, to=1, length=100)
coef <- sapply(p.shared, tiling.shared)
plot(p.shared, coef, pch=20, main='This should monotonically increase')
expect_true( abs(coef[1]) < 0.1)
expect_equal( coef[length(coef)], 1)



context("Anti-correlated trains.")
tiling.anti <- function(dt=0.05, n=2000) {
  master <- seq(from=0, by=0.5, length=n)
  odd <- rep_len(c(TRUE, FALSE), n)
  a <- master[odd]
  b <- master[!odd]
  tiling.corr(a,b,dt)
}

dts <- seq(from=0.05, to=1.0, length=100)
dts <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
coeff <- sapply(dts, tiling.anti)
plot(dts, coeff, pch=20, type='b')
## We get NaN for dt=0.5 and dt=1 second.




context("Checking the rec.time works")

## Generate a pair of uncorrelated trains, a and b.
## Then make trains a' and b' by simply adding a constant Z to all times.
## then check that tiling(a,b) == tiling(a+z, b+z)
beg <- 0; end <- 2000; n <- 3000;
z <- 5000;                              #large offset for 2nd set of trains
for (i in seq_len(10)) {
  a <- sort(runif(n, beg, end))
  b <- sort(runif(n, beg, end))
  c1 <- tiling.corr(a, b, rec.time=c(beg, end))
  c2 <- tiling.corr(a+z, b+z, rec.time=z + c(beg, end))
  all.equal(c1, c2)
}

context("Pathological corner case with synthetic trains")
## This is when Pa=Tb=1 so both numerator and denominator are zero.
## What should we do about this case?  Unlikely to happen for
## realistic trains.
a <- 1; b <- 2                          # one spike in each time.
expect_equal(tiling.corr(a, b, dt=1, rec.time=c(0, 3)), 1)
expect_equal(tiling.corr(a, b, dt=1), NaN) #is this correct?!?


context("Check the array wide computation")


data.file <- system.file("examples", "P9_CTRL_MY1_1A.txt",
                         package = "sjemea")
s <- jay.read.spikes(data.file)
system.time(t1 <- tiling.allpairwise.old(s))
system.time(t2 <- tiling.allpairwise(s))
require(lattice)
levelplot(t1)
levelplot(t2)
u <- upper.tri(t1, diag=TRUE)
expect_equal(t1[u], t2[u])

