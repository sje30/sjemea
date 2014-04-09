
tiling.allpairwise <- function(s, dt=0.05) {
  ## Return matrix of all pairwise tiling correlations from object s.
  ## TODO: only need upper triangular elements completed.
  N <- s$NCells
  indices=array(0,dim=c(N,N))

  for(i in 1:N) {
    for(j in 1:N) {

      z <- .C("run_TM",as.integer(length(s$spikes[[i]])),
              as.integer(length(s$spikes[[j]])),
              as.double(dt),
              as.double(s$rec.time[[2]]),
              index=as.double(1),
              as.double(as.vector(s$spikes[[i]])),
              as.double(as.vector(s$spikes[[j]])),
              PACKAGE="sjemea")
      indices[i,j] <- z[[5]]
    }
  }
  indices
}

