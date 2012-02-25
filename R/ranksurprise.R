## ranksurprise.R --- Rank Surprise method for burst detection.
## Author: Zhengzheng Zhang
## Copyright: GPL.
## Following method in Gourevitch and Eggermont (2007)

## "Rank Surprise" method for burst detection. The method and algorithm are 
## described in the paper, a nonparametric approach for detection of bursts 
## in spike trains. See reference.


val2rk <- function(values) {
  ## Convert values to ranks, with mean of ranks for tied values. Alternative
  ## and faster version of "tiedrank" in statistical toolbox of Matlab 6-7.
  
  lp = length(values)
  rk = rk2 = rep(0,lp)
  S = sort(values, index.return = TRUE)
  y = S$x
  cl = S$ix
  rk[cl] = c(1:lp)
  cl2 = sort(-values, index.return = TRUE)$ix
  rk2[cl2] = c(1:lp)
  ranks = (lp+1-rk2+rk)/2
  ranks
}


rs.find.bursts <- function(s, limit = NULL, RSalpha = -log(0.01)) {
  nspikes = length(s)
  if (nspikes <= 10) {
    archive_burst = NA
  }else{
    ## General parameters
    
    ## limit for using the real distribution
    q_lim = 30
    ## minimum length of a burst in spikes
    l_min = 3
    
    ## General vectors

    ## vector (-1)^k
    alternate = rep(1, 400)
    alternate[c(1:200)*2] = -1
    ## log factorials
    log_fac = cumsum(log(1:q_lim))
    
    ## Ranks computation
    
    ## compute the ISI
    ISI = diff(s)
    N = length(ISI)
    ## ISI value not to include in a burst
    if (length(limit)==0){
      ## percentile 75% (default)
      limit=summary(ISI)["3rd Qu."][[1]]
    }

    ## compute ranks
    R = val2rk(ISI)
    
    ## Find sequences of ISI under 'limit'
    
    D = rep(0, N)
    D[ISI<limit]=1
    ISI_limit = diff(D)
    ## first time stamp of these intervals
    begin_int = which(ISI_limit==1)+1
    ## manage the first ISI
    if (ISI[1]<limit){
      begin_int = c(1, begin_int)      # the first IS is under limit
    }    
    ## last time stamp of these intervals
    end_int = which(ISI_limit==-1)
    ## manage the last ISI
    if (length(end_int)<length(begin_int)){
      end_int=c(end_int, N)
    }
    ## length of intervals of interest
    length_int=end_int-begin_int+1
    
    ## Initializations
    archive_burst_RS=c()
    archive_burst_length=c()
    archive_burst_start=c()
    archive_burst_end=c()
    archive_burst_IBI=c()
    archive_burst_durn=c()
    archive_burst_mean.isis=c()

    ## Going through the intervals of interest
    indic=0
    for (n_j in begin_int){
      indic=indic+1
      p_j=length_int[indic]
      subseq_RS = NULL 
      subseq_RS_RS=c()
      subseq_RS_i =c()
      subseq_RS_q =c()     
      ## test each set of spikes
      for (i in 0:(p_j-(l_min-1))){
        ## length of burst tested
        q=l_min-2
        while (q<(p_j-i)){
          q=q+1
          ## statistic
          u=sum(R[(n_j+i):(n_j+i+q-1)])
          u=floor(u)
          if (q<q_lim){
            ## exact discrete distribution
            k=c(0:((u-q)/N))
            length_k=length(k) 
            SUM = 0                     
            for (j in 1:q){
              SUM = SUM + log(u-matrix(rep(k,q),q,length_k, byrow = TRUE)*N
                - matrix(rep(c(0:(q-1)),length_k),length(c(0:(q-1))),length_k))[j,]
            }      
            if (length_k<2){
              prob = exp((SUM- log_fac[1]-log_fac[q-k])-q*log(N))%*%alternate[1:length_k]
            }else{
              prob = exp((SUM- log_fac[c(1, k[2:length_k])]-log_fac[q-k])-q*log(N))%*%alternate[1:length_k]
            }    
          }else{
            ## approximate Gaussian distribution
            prob=pnorm((u-q*(N+1)/2)/sqrt(q*(N^2-1)/12))
          }
          RS=-log(prob)
          ## archive results for each subsequence [RSstatistic beginning length]
          if (RS>RSalpha){
            subseq_RS_RS[length(subseq_RS_RS)+1]=RS
            subseq_RS_i[length(subseq_RS_i)+1]=i
            subseq_RS_q[length(subseq_RS_q)+1]=q
            subseq_RS = rbind(subseq_RS_RS, subseq_RS_i, subseq_RS_q)
            subseq_RS = t(as.matrix(subseq_RS))
          }
        }
      }
      ## vet results archive to extract most significant bursts
      if (length(subseq_RS)!=0){
        ## sort RS for all subsequences
        if (!is.vector(subseq_RS)){
          ind = sort(subseq_RS[,1], decreasing = TRUE, index = TRUE)$ix
          subseq_RS = subseq_RS[ind,]
        }           
        while (length(subseq_RS)!=0){
          ## extract most surprising burst
          if (is.vector(subseq_RS)){
            current_burst=subseq_RS
            archive_burst_RS[length(archive_burst_RS)+1]=current_burst[1]
            archive_burst_length[length(archive_burst_length)+1]=current_burst[3]+1  #number of ISI involved + 1
            archive_burst_start[length(archive_burst_start)+1]=n_j+current_burst[2]
            archive_burst_end[length(archive_burst_end)+1]=n_j+current_burst[2]+current_burst[3]
            subseq_RS=NULL
          }else{                      
            current_burst=subseq_RS[1,]
            archive_burst_RS[length(archive_burst_RS)+1]=current_burst[1]
            archive_burst_length[length(archive_burst_length)+1]=current_burst[3]+1 #number of ISI involved + 1
            archive_burst_start[length(archive_burst_start)+1]=n_j+current_burst[2]
            archive_burst_end[length(archive_burst_end)+1]=n_j+current_burst[2]+current_burst[3]
            ## remove most surprising burst from the set
            ## subseq_RS=subseq_RS(2:end,:);
            ## keep only other bursts non-overlapping with this burst
            D = (subseq_RS[,2]+subseq_RS[,3]-1)<current_burst[2]
            E = subseq_RS[,2]>(current_burst[2]+current_burst[3]-1)
            subseq_RS=subseq_RS[D|E,]
          }
        }
      }
    }
    
    ## sort bursts by ascending time
    if (length(archive_burst_RS)!=0){
      S = sort(archive_burst_start, index.return = TRUE)
      archive_burst_start = S$x
      ind_sort = S$ix
      archive_burst_RS=archive_burst_RS[ind_sort]
      archive_burst_length=archive_burst_length[ind_sort]
      archive_burst_end=archive_burst_end[ind_sort]
      last.end = archive_burst_end[-length(archive_burst_end)]
      last.end = c(NA, s[last.end])
      archive_burst_IBI = s[archive_burst_start]-last.end
      archive_burst_durn = s[archive_burst_end]-s[archive_burst_start]
      archive_burst_mean.isis = archive_burst_durn/(archive_burst_length-1)
      archive_burst=cbind(beg=archive_burst_start, end=archive_burst_end, IBI=archive_burst_IBI,
        len=archive_burst_length, durn=archive_burst_durn, 
        mean.isis=archive_burst_mean.isis, SI=archive_burst_RS)
    }else{
      archive_burst = NA
    }
  }
  archive_burst
} 
