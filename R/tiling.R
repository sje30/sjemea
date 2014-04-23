
tiling.allpairwise <- function(s, dt=0.05) {
  ## Return matrix of all pairwise tiling correlations from object s.
  ## TODO: only need upper triangular elements completed.
  N <- s$NCells
  indices=array(0,dim=c(N,N))

  duration <- as.double( diff(s$rec.time))
  for(i in 1:N) {
    ni <- s$nspikes[i]
    for(j in 1:N) {
      nj <- s$nspikes[j]
      z <- .C("run_TM",
              as.integer(ni),
              as.integer(nj),
              as.double(dt),
              duration,
              index=double(1),
              as.double(s$spikes[[i]]),
              as.double(s$spikes[[j]]),
              PACKAGE="sjemea")
      indices[i,j] <- z$index
    }
  }
  indices
}

