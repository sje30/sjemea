#include <R.h>
#include <S.h>			/* for seed_in, seed_out */
#include <Rinternals.h>
#include <R_ext/Rdynload.h>	/* for DLL types */


void run_TM(int *N1v,int *N2v,double *dtv,double *Tv,double *index,double *spike_times_1,double *spike_times_2);

void ns_count_activity(Sfloat *allspikes, int *nspikes, int *pncells,
		       Sfloat *pbeg, Sfloat *pend, Sfloat *pwid,
		       int *pnbins,
		       int *count);

void frate(Sfloat *allspikes, int *nspikes, int *pncells,
	   Sfloat *pbeg, Sfloat *pend, Sfloat *pwid,
	   int *pnbins,
	   double *counts);

void count_overlap_arr(Sfloat *spikes,
		       int *pn,
		       int *nspikes,
		       int *first_spike,
		       int *rates_ok,
		       int    *pno_min,
		       Sfloat *pt, /* duration of recording */
		       Sfloat *pdt,
		       Sfloat *corrs /* return array */);

void coincident_arr(Sfloat *a, int *pna,
		    Sfloat *bs, int *nb, int *pnchannels,
		    int *close, Sfloat *pw);

void bin_overlap(Sfloat *a, int *pna, Sfloat *b, int *pnb, Sfloat *pdt,
		 int *bins, int *pnbins);

void bin2_overlap(Sfloat *a, int *pna, Sfloat *b, int *pnb, Sfloat *pdt,
		  int *bins, int *pnbins);

void count_overlap(Sfloat *a, int *pna, Sfloat *b, int *pnb, Sfloat *pdt,
		   int *res);

void tiling_arr(Sfloat *spikes,
		int *pn,
		int *nspikes,
		int *first_spike,
		int *rates_ok,
		int    *pno_min,
		Sfloat *rec_time, /* recording time */
		Sfloat *pdt,
		Sfloat *corrs /* return array */);


static R_NativePrimitiveArgType run_TM_t[7] =
  {INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,  REALSXP, REALSXP};

static R_NativePrimitiveArgType ns_count_activity_t[8] =
  {REALSXP, INTSXP, INTSXP,
   REALSXP, REALSXP, REALSXP,
   INTSXP, INTSXP};

static R_NativePrimitiveArgType frate_t[8] =
  {REALSXP, INTSXP, INTSXP,
   REALSXP, REALSXP, REALSXP,
   INTSXP, REALSXP};


static R_NativePrimitiveArgType count_overlap_arr_t[9] =
  {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
   REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType bin_overlap_t[7] =
  {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP,
   INTSXP, INTSXP};

static R_NativePrimitiveArgType bin2_overlap_t[7] =
  {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP,
   INTSXP, INTSXP};

static R_NativePrimitiveArgType coincident_arr_t[7] =
  {REALSXP, INTSXP,
   REALSXP, INTSXP, INTSXP,
   INTSXP, REALSXP};

static R_NativePrimitiveArgType count_overlap_t[6] =
  {REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType tiling_arr_t[7] =
  {REALSXP, INTSXP, INTSXP, INTSXP, 
   REALSXP, REALSXP, REALSXP};

/* Let's register the functions here. */
R_CMethodDef cMethods[] = {
  {"run_TM",            (DL_FUNC) &run_TM, 7, run_TM_t},
  {"ns_count_activity", (DL_FUNC) &ns_count_activity, 8, ns_count_activity_t},
  {"frate",             (DL_FUNC) &frate, 8, frate_t},
  {"count_overlap_arr", (DL_FUNC) &count_overlap_arr, 9, count_overlap_arr_t},
  {"coincident_arr",    (DL_FUNC) &coincident_arr, 7, coincident_arr_t},
  {"bin_overlap",       (DL_FUNC) &bin_overlap, 7, bin_overlap_t},
  {"bin2_overlap",      (DL_FUNC) &bin2_overlap, 7, bin2_overlap_t},
  {"count_overlap",     (DL_FUNC) &count_overlap, 6, count_overlap_t},
  {"tiling_arr",        (DL_FUNC) &tiling_arr, 7, tiling_arr_t},
  {NULL, NULL, 0}
};


void R_init_sjemea(DllInfo *info) {
  /* Register the routines that you want to call from R. */
  /* This should me name R_init_XXX whwere XXX is your package name. */
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}


   

/* We use this SMALLVAL to test whether two floating point values can
 * be regarded as "equal": we see if the absolute difference between
 * the two values is less than this SMALLVAL. */
#define SMALLVAL 1e-12


/* Theory of binning numbers.  Seems quite straightforward, but there
 * is a little thing to worry about it: numbers falling onto a bin
 * edge.  If we want to bin a set of numbers that vary between MIN and
 * MAX into n bins, we make each bin of width w:
 *   w = (MAX-MIN) / n
 * and then for a value x, we increment bin b, where
 *   b =  (int) ( (x - MIN)/ w)
 * where (int) rounds down to nearest integer.  This makes each bin of
 * the type [low, high).  So, for MIN=0, MAX=3, n=3, the overall range
 * is divided into say 3 bins like: [0, 1) [1, 2) [2, 3).  This means
 * that the overall range of the histogram will be [0,3).  So, if we
 * try to bin a value of x=3 (the max value), it falls outside the
 * last bin.  In R, this problem is overcome using the INCLUDE.LOWEST
 * variable.  In this C code, we explicitly test for this case, so the
 * max value is included in the last bin.
 *
 * For more information on binning, see David Young's help for the
 * POP-11 routine array_hist, included at the end of this file.
 */
  
void count_overlap(Sfloat *a, int *pna, Sfloat *b, int *pnb, Sfloat *pdt,
		  int *res)
{
  /* A[i] is the time of the ith spike in cell A.  A has *PNA spikes.
   * Likewise for cell B.  Assume that spike times are ordered,
   * earliest first.  We return (in *RES) the number of spikes in B
   * that occur within +/- dt of a spike in A.
   */

  int count=0;
  int sa, sb, low;
  Sfloat alow, ahigh, dt;
  int na, nb;

  na = *pna; nb = *pnb; dt=*pdt;

  /* low is the index of first spike time in train b that could be
   * withing the range [a-dt, a+dt] of the spike in a.  We rely on the
   * spike times being ordered to use this trick.
   */
   
  low = 0;
  for( sa=0; sa < na; sa++) {
    alow = a[sa]-dt; ahigh= a[sa]+dt;
    for(sb=low; sb < nb;sb++) {
      if (b[sb] > ahigh)
	sb = nb;		/* stop checking this train */
      else			/* spike could be in [alow, ahigh] */
	if (b[sb] >= alow)	/* equality holds when delta_t is -tmax */
	  count++;
        else
	  low = sb;		/* can ignore this spike in B from now on. */
    }
  }
  *res = count;
}

void count_overlap_arr(Sfloat *spikes,
		       int *pn,
		       int *nspikes,
		       int *first_spike,
		       int *rates_ok,
		       int    *pno_min,
		       Sfloat *pt, /* duration of recording */
		       Sfloat *pdt,
		       Sfloat *corrs /* return array */) {

  int a, b, n, no_min, count;
  Sfloat *sa, *sb; 		/* pointers to current spike trains  */
  int n1, n2;
  Sfloat k1, res;
  int debug;

  debug = 0;
  
  n = *pn; no_min=*pno_min;
  res = -999;
  k1 = *pt / (2.0 * *pdt);		/* simple constant for correlation index. */
  
  if (debug) {
    Rprintf("Time %.2f Deltat %.3f no_min %d k1 %.2f\n",
	    *pt, *pdt, no_min, k1);
  }

	    
  for (a=0; a<n-1; a++) {
    n1 = nspikes[a];
    sa = &(spikes[first_spike[a]]);
    
    for (b=a+1; b<n; b++) {
      n2 = nspikes[b];
      sb = &(spikes[first_spike[b]]);
      
      if (no_min || (rates_ok[a] && rates_ok[b])) {
	if (debug)
	  Rprintf("first spikes %.2f %.2f n %d %d\n",
		  *sa, *sb, n1, n2);
	
	count_overlap(sa, &n1, sb, &n2, pdt, &count);
	/* res = (count*k1) / ((double)(n1 * n2 ));*/
	res = (count*k1) / ((double)(n1) * (double)(n2));
	if (debug)
	  Rprintf("%d %d: count %d corr %.2f\n", a, b, count, res);
      } else {
	/* One of the spikes below min firing rate. */
	res = R_NaN;		/* defined in R_ext/Arith.h */
      }

      corrs[(b*n)+a] = res;
      
    }
  }
}

void bin_overlap(Sfloat *a, int *pna, Sfloat *b, int *pnb, Sfloat *pdt,
		 int *bins, int *pnbins)
{
  /* bin_overlap is similar to count_overlap except that we return a
   * histogram which bins the time difference between spikes into one
   * of several time bins.  The maximum absolute time difference is
   * *PDT; the histogram *INS returned is on length *PBINS.  Here we
   * ignore whether the time difference is +ve or -ve.  To check that
   * the binning is working, run the R function test.count.hist.nab().
   * Each histogram bin is of the form [low1, high1) with the last bin
   * specially set to [lown,highn], so the overall range of this
   * histogram is [0,T], where T=*PDT.
   */
  
  int sa, sb, low;
  Sfloat alow, ahigh, dt, delta_t, bin_wid, max_val;
  int na, nb;
  int bin_num, nbins;

  na = *pna; nb = *pnb; dt=*pdt; nbins = *pnbins;
  max_val = dt;
  bin_wid = dt/(double)nbins;
  /* low is the index of first spike time in train b that could be
   * withing the range [a-dt, a+dt] of the spike in a.  We rely on the
   * spike times being ordered to use this trick.
   */
   
  low = 0;
  for( sa=0; sa < na; sa++) {
    alow = a[sa]-dt; ahigh= a[sa]+dt;
    for(sb=low; sb < nb;sb++) {
      if (b[sb] > ahigh)
	sb = nb;		/* stop checking this train */
      else			/* spike could be in [alow, ahigh] */
	if (b[sb] >= alow) {	/* equality when (deltat == -tmax) */
	  /* need to bin this value. */
	  delta_t = fabs( b[sb] - a[sa]);
	  bin_num = (int) (delta_t / bin_wid);
	  if ((bin_num == nbins) &&(fabs(delta_t - max_val) < SMALLVAL))
	    bin_num--;		/* put max value into largest bin. */
	  if ( (bin_num <0 ) || (bin_num >= nbins))
	    Rprintf("bin number wrong %f %d\n", delta_t, bin_num);
	  else {
	    bins[bin_num]++;
	  }
	} else
	  low = sb;		/* can ignore this spike in B from now on. */
    }
  }
}


void bin2_overlap(Sfloat *a, int *pna, Sfloat *b, int *pnb, Sfloat *pdt,
		 int *bins, int *pnbins)
{
  /* bin2_overlap is a bidirectional version of bin_overlap.
   * This time the sign of the time difference between spikes is important.
   * Each histogram bin is of the form [low1, high1) with the last bin
   * specially set to [lown,highn], so the overall range of this
   * histogram is [-T,T] where T=*PDT.
   */

  int sa, sb, low;
  Sfloat alow, ahigh, dt, delta_t, bin_wid, min_val, max_val;
  int na, nb;
  int bin_num, nbins, bin_numi;

  na = *pna; nb = *pnb; dt=*pdt; nbins = *pnbins;
  min_val = 0.0-dt;		/* smallest value that we will bin. */
  max_val = dt;
  /* the range of times is now [-dt,dt], so we need to double dt. */
  bin_wid = (2.0 * dt)/nbins;
  /*Rprintf("bin_wid %f\n", bin_wid);*/
  /* low is the index of first spike time in train b that could be
   * withing the range [a-dt, a+dt] of the spike in a.  We rely on the
   * spike times being ordered to use this trick.
   */
   
  low = 0;
  for( sa=0; sa < na; sa++) {
    alow = a[sa]-dt; ahigh= a[sa]+dt;
    for(sb=low; sb < nb;sb++) {
      if (b[sb] > ahigh)	/* >= for case when at max. */
	sb = nb;		/* stop checking this train */
      else			/* spike could be in [alow, ahigh] */
	if (b[sb] >= alow) {	/* equality when (deltat == -tmax) */
	  /* need to bin this value. */
	  delta_t = ( b[sb] - a[sa]);

	  /* Compute bin number using both floor and casting to an
	   * int.  This is temporary code as I think I found that
	   * sometimes the "self spike" for auto-correlations was put
	   * in the wrong bin.  If we get this error, we are in
	   * trouble.
	   */
	  bin_numi = (int)((delta_t - min_val)/ bin_wid);
	  bin_num = (int)floor((delta_t - min_val)/ bin_wid);
	  /*
	  if (bin_num != bin_numi) {
	    Rprintf("XXX different bin numbers: dt %f min %f wid %f floor %d (int) %d %f\n",
		    delta_t, min_val, bin_wid,
		    bin_num, bin_numi, delta_t);
	  }
	  */
	  if ( (bin_num == nbins) &&(fabs(delta_t - max_val) < SMALLVAL))
	    bin_num--;		/* fits into largest bin. */
	  if ( (bin_num <0 ) || (bin_num >= nbins))
	    Rprintf("bin2: number wrong %f %d %f\n",delta_t, bin_num, bin_wid);
	  else {
	    /* When making auto-correlation, see which bin the
	     * "self-comparison" is put into:
	     */
	    /* if (sa==sb) */
	    /*Rprintf("%d: when binning own spike (time %f), dt %f bin %d\n",*/
	    /*sa, a[sa], delta_t, bin_num);*/

	    bins[bin_num]++;
	  }
	} else
	  low = sb;		/* can ignore this spike in B from now on. */
    }
  }
}


void ns_count_activity(Sfloat *allspikes, int *nspikes, int *pncells,
		       Sfloat *pbeg, Sfloat *pend, Sfloat *pwid,
		       int *pnbins,
		       int *count)
{

  /* Compute the network spike activity.
   *
   * ALLSPIKES = vector of spike times, flattened, so that we get all
   * spikes for unit 1, then unit 2, and so on
   *
   * NSPIKES[j] indicates the number of spikes from cell j.
   *
   * NCELLS = number of cells.
   * BEG, END = first and last spike time.
   * WID = duration of each network spike bin.
   * NBINS = number of bins (calculated in R, rather than C)
   *
   * COUNT[i] = stores the number of units that were active in bin i.
   * If a spike train fires more than once during a bin, it counts
   * only once (using the LAST flag below).
   */
  
  Sfloat *p, beg, end, wid;
  int ncells, last, b, n, unit, nbins;
  
  ncells = *pncells; beg = *pbeg; end = *pend; wid = *pwid;
  nbins = *pnbins;

  p = allspikes;
  for (unit=0; unit<ncells; unit++) {
    /* Count the spikes on electrode UNIT. */
    n = nspikes[unit];
    last = -1;			/* check to only increment bin once per unit. */

    while(n-- >0) {
      b = (int) ( (*p++ - beg)/wid); /* calc bin number; increment spike ptr */

      /* Check bin number is valid: shouldn't happen. */
      if ( (b <0 ) || (b >= nbins))
	Rprintf("bin number wrong %f %d\n", *(p-1), b);
      else {
	/* Update count in relevant bin. */
	if (last != b) {
	  count[b]++;
	  last = b;		/* stop this bin being updated again for
				 * current unit. */
	}
      }
    }
  }

}

void frate(Sfloat *allspikes, int *nspikes, int *pncells,
	   Sfloat *pbeg, Sfloat *pend, Sfloat *pwid,
	   int *pnbins,
	   double *counts)
{

  /* Compute the firing rate.
   * Adapted from ns_count_activity.
   *
   * ALLSPIKES = vector of spike times, flattened, so that we get all
   * spikes for unit 1, then unit 2, and so on
   *
   * NSPIKES[j] indicates the number of spikes from cell j.
   *
   * NCELLS = number of cells.
   * BEG, END = first and last spike time.
   * WID = duration of each network spike bin.
   * NBINS = number of bins (calculated in R, rather than C)
   *
   * COUNTS[i,c] = for cell c, estimate the firing rate of the ITH bin.
   * COUNTS is a long vector of size (NBINS*NCELLS)
   */
  
  double *p, beg, end, wid, *count;
  int ncells, b, n, unit, nbins, skip;
  
  ncells = *pncells; beg = *pbeg; end = *pend; wid = *pwid;
  nbins = *pnbins;

  p = allspikes;
  count = counts;
  skip =0;
  
  for (unit=0; unit<ncells; unit++) {
    /* Count the spikes on electrode UNIT. */
    n = nspikes[unit];

    while(n-- >0) {
      b = (int) ( (*p++ - beg)/wid); /* calc bin number; increment spike ptr */

      /* Check bin number is valid: shouldn't happen. */
      if ( (b <0 ) || (b >= nbins))
	/* Rprintf("bin number wrong %f %d\n", *(p-1), b); */
	skip++;
      else {
	/* Update count in relevant bin. */
	  count[b]++;
      }
    }
    /* finished checking one electrode, so prepare for next electrode.*/
    count += nbins;
  }

  /* After processing all electrodes, divide the counts by the bin
   * width to estimate firing rate.
   * TODO
   */
}


void arraywide_autocorr(Sfloat *allspikes, int *nspikes, int *pncells,
			Sfloat *pwid,
			int *pnbins,
			int *count)
{

  /* Compute the network spike activity.
   *
   * ALLSPIKES = vector of spike times, flattened, so that we get all
   * spikes for unit 1, then unit 2, and so on
   *
   * NSPIKES[j] indicates the number of spikes from cell j.
   *
   * NCELLS = number of cells.
   * WID = duration of each autocorrelation bin.
   * NBINS = number of bins.
   *
   * COUNT[i] = stores the number of times that two spikes on the same
   * train were within a certain time apart (of width WID).
   */
  
  Sfloat wid, s_i;
  int ncells, b, n, unit, nbins, i, j, looking;
  int first_spike, last_spike;		/* first, last spike index of any train. */
  ncells = *pncells; wid = *pwid;
  nbins = *pnbins;

  first_spike = 0;

  for (unit=0; unit<ncells; unit++) {
    /* Compute autocorrelation  for spikes on electrode UNIT. */
    n = nspikes[unit];

    last_spike = first_spike + n;
    
    for (i=first_spike; i<last_spike-1; i++) {
      /* Autocorrelate spike i with "future" spikes.  No need to check
       * let i equal the last spike of a train, since there are no
       * future spikes on the train to correlated with. */
      s_i = allspikes[i];
      j = i+1;
      looking = TRUE;
      while (looking) {
	b = (int) ( ( allspikes[j] - s_i)/wid);
	/*assert(b>=0);*/
	if (b > nbins) {
	  looking = FALSE;
	} else {
	  count[b]++;
	  j++;
	  if (j==last_spike) {
	    looking = FALSE;
	  }
	}
      }
    }

    /* Finished checking all spikes on current train; move onto next train. */
    first_spike += n;
  }

}



void coincident(Sfloat *a, int *pna, Sfloat *b, int *pnb, int *close,
		 Sfloat *pw) {
  int i, j;
  int spikes_in_b, looking, res;
  int n_a, n_b;
  Sfloat diff, t_a, w;

  w = *pw;
  i = 0; j = 0;
  n_a = *pna; n_b = *pnb;
  
  spikes_in_b = (n_b > 0);
  while ( (i < n_a) && spikes_in_b) {

    t_a = a[i];
    looking = TRUE; res=0;
    while (looking) {
      diff = b[j] - t_a;

      if ( diff < -w) {
	/* case 1 - spike j too early for any event in A */
	j++;
	if (j == n_b) {
	  /* no more spikes in b, so end everything */
	  spikes_in_b = FALSE; looking = FALSE;
	}
      } else {
	if ( diff < w) {
	  /* case 2 - good, coincident */
	  res = 1; looking = FALSE;
	} else {
	  /* must be case 3 - spike j too late for current event in A */
	  /* when debugging, can set res=-1 to see this case */
	  res = 0; looking = FALSE;
	}
      }

    }
    /* end checking for time i in A */
    close[i] = res;
    i++;
  }
}

void coincident_arr(Sfloat *a, int *pna,
		    Sfloat *bs, int *nb, int *pnchannels,
		    int *close, Sfloat *pw)
{
  /* Compute overlap across all channels of the spikes.
   * It calls "coincident" for each spike train in BS.
   * 
   * A = vector (of length *PNA) of all reference times, e.g. times
   * of peak network spike activity.
   * BS = all  spikes from all channels collapsed into one long
   * vector (first, all the spikes for channel 1, then all spikes for
   * channel 2, and so on.)
   * NCHANNELS = number of channels.
   * NB[j] = number of spikes in channel J.  The sum of elements in NB should
   * equal the total number of spikes, i.e, the length of BS.
   * CLOSE = long vector of output, coerced into an array by R.
   * PW = window size.
   */
  int nchann;
  int j, n_a, n_b;
  Sfloat *b;


  /* We use pointer arithmetic in two ways here:
   * 1. B is incrementally updated to move to the start of the next
   * spike train within BS.
   * 2. CLOSE is a vector of length (N_A * N_CHANNELS); for each
   * channel we fill in N_A values, and so for each channel, CLOSE is
   * moved on by N_A elements to gradually fill up the whole vector.
   * R then converts the vector into a matrix. */
   
  n_a = *pna;
  b = bs;			/* point to start of current spike train */
  nchann = *pnchannels;
  for (j = 0; j < nchann; j++) {
    n_b = nb[j];		/* #spikes in current train */
    coincident(a, pna, b, &n_b, close, pw);
    close += n_a;		/* move to next free bit of output array */
    b += n_b;			/* move onto next spike train */
  }



}

#define MAX(X, Y)  ((X) > (Y) ? (X) : (Y))
#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

void bci_calc(int *pn, Sfloat *beg, Sfloat *end,
	      int *nbursts, int *start,
	      Sfloat *durn,
	      Sfloat *res) {
  int n;
  int a, b;

  Sfloat olap, this_olap;
  Sfloat beg_a, end_a, beg_b, end_b;
  int ind_a, ind_b, a1, b1;
  int debug=0;
  
  n = *pn;			/* number of cells */
  Rprintf("%d cells\n", n);
  for (a = 0; a <n-1; ++a) {
    for (b = a+1; b < n; ++b) {
      if (debug) Rprintf("%d %d\n", a, b);

      /* Compute overlap between cells A and B. */
      olap=0.0;
      ind_a = start[a];
      for (a1=nbursts[a]; a1>0; a1--) {
	beg_a = beg[ind_a];
	end_a = end[ind_a];

	/* Could also optimize not to start at start of bursts for b,
	 * but use information from previous comparisons as to where to start.
	 */

	/* Look over bursts in b. */
	ind_b = start[b];
	for (b1=nbursts[b]; b1>0; b1--) {
	  beg_b = beg[ind_b];
	  end_b = end[ind_b];

	  if (beg_b > end_a) {
	    /* No more bursts in b can overlap with current burst in a. */
	    b1=0;
	  } else {
	    this_olap = MIN(end_a, end_b) - MAX(beg_a, beg_b);
	    if (debug) {
	      Rprintf("%.2f %.2f  vs %.2f %.2f => %.2f\n",
		      beg_a, end_a, beg_b, end_b, this_olap);
	    }
	    olap += MAX(this_olap, 0);
	  }
	  ind_b++;
	}
	ind_a++;
      }
      /* Store results. */
      res[a + (b*n)] = olap / durn[b]; /* b; */
      res[b + (a*n)] = olap / durn[a]; /* 0 - a; */
    }
  }
}




/**********************************************************************/
/*
HELP ARRAY_HIST                             David Young, January 1994


LIB *array_hist provides a procedure for obtaining a histogram of the
values in a region of an array.  All the values must be numbers.

         CONTENTS - (Use <ENTER> g to access required sections)

 -- Procedure array_hist
 -- Counting integer values
 -- Counting floating point values
 -- External optimisation
 -- Re-using histogram vectors
 -- Offsetting results in the vector

-- Procedure array_hist -----------------------------------------------

array_hist(__array, region, low, __nbins, _high) -> (_nlow, _hist, __nhigh)

        The histogram is formed for values in the part of __array
        specified by region.  The list region is in *boundslist style
        (i.e. first two elements give range of indices in first
        dimension, next two give range in second dimension, etc.). If
        region is false, the whole of the array is examined.

        The numbers low and _high give the overall range of values to
        count in the histogram. (The procedure *array_mxmn may be useful
        for obtaining these in the general case.) The range between low
        and _high is divided into __nbins equal parts.

        The bin width (the range of values that get counted in one bin)
        is given by

            __binwidth = (_high - low) / __nbins

        The result _hist is a vector containing counts of the values in
        each bin. The results _nlow and __nhigh return the counts of values
        that fell outside the range covered by _hist.

        To be precise, low is the smallest value to get counted in the
        first bin and _high is the smallest value _just too __large to get
        counted in the last bin. (This means that the treatment of
        integers and floats can be consistent.) A value _V from the array
        is treated as follows:

            _V < low:        increment _nlow
            _V >= _high:      increment __nhigh
            otherwise:      increment _hist(_I) where

                 _I = floor( (_V - low) * __nbins / (high - low) ) + 1

        Apart from rounding errors, this means that in the last case _I
        is chosen such that

                low + __binwidth * (_I - 1) <= _V < low + __binwidth * _I

        (The floor function returns the largest integer less than or
        equal to its argument.)

        It is possible to re-use vectors and to place the counts in the
        vector starting from some element other than the first. These
        options are described below.

-- Counting integer values --------------------------------------------

Suppose the values in __array are integers in the range 0 ... 255, and we
want to know how many of each there are.  The correct call is

    array_hist(array, false, 0, 256, 256) -> (nlow, hist, nhigh);

The __nbins argument is 256 because there are 256 different values to
count. Note that _high is 256, not 255, because it must be the next value
above the top of the histogram range. To make the bin width equal to 1,
we need

    _high = low + __nbins

The element _hist(_I) will contain the number of values in the array equal
to _I-1.  The -1 is necessary because the values start at 0, but vectors
are indexed from 1.

In general, for integer values to be counted properly, with _K different
values counted in each bin, we need

    _high = low + _K * __nbins

and the _I'th element of _hist will contain the count for values in the
range low + _K * (_I-1) to low + (_K+1) * (_I-1) - 1.

To sum up, to count integers in the range __N0 to __N1 inclusive you should
use:

    low = __N0
    _high = __N1 + _1

and the number of different values counted in each bin will be

    _K = (_high - low) / __nbins

with __nbins chosen to make _K an integer.

For example, we can look at the performance of the POP-11 random number
generator by filling an array with random numbers in the range 1 to 16
and looking at its histogram.

    vars arr, nlo, hist, nhi;
    newarray([1 1000], erase <> random(% 16 %)) -> arr;
    array_hist(arr, false, 1, 16, 17) -> (nlo, hist, nhi);

    nlo =>
    ** 0
    hist =>
    ** {59 62 67 75 60 59 66 61 56 67 47 64 60 70 58 69}
    nhi =>
    ** 0

As expected, no values are less then 1 or greater than or equal to 17,
and the 1000 values are reasonably evenly distributed. (You will not get
an identical distribution if you try this.)

-- Counting floating point values -------------------------------------

Counting floating point values ("decimals" in POP-11) is usually
simpler, as low and _high then normally correspond exactly to the range
of interest.  For example, to test the performance of the random number
generator on floats, we can fill an array with numbers from 0.0 to 1.0
and look at its histogram in much the same way as before:

    newarray([1 1000], erase <> random0(% 1.0 %)) -> arr;
    array_hist(arr, false, 0.0, 16, 1.0) -> (nlo, hist, nhi);

The results will be similar to the previous example. The bin width in
this case is 0.0625 - sixteen of these cover the range from 0.0 to 1.0.

Rounding errors mean that values on, or very close to, bin boundaries
may get counted in the wrong bin.  This risk is inevitable with floating
point calculations.  If the values fall into natural groups, the problem
can be eliminated by putting the bin boundaries firmly into the gaps.
For example, if the values are whole numbers (although represented as
floats) in the range A0 to A1 inclusive, and the bin width is to be 1,
then it would be sensible to use

    low = __A0 - 0.5
    _high = __A1 + 0.5
    __nbins = round(_high - low)

However, this should not be done if the values are actually represented
as integers - see the section above.

-- External optimisation ----------------------------------------------

Two cases are dealt with using external code, for much increased speed:

    1. __array is a packed array of single precision floating point
    values, as produced for example by *newsfloatarray.

    2. __array is a packed array of bytes, as produced for example by
    *newbytearray, and both low and _high are integers.

The result _hist will be an *INTVEC.

-- Re-using histogram vectors -----------------------------------------

It is possible to re-use a histogram vector to avoid creating garbage,
by passing it as an argument as __nbins (instead of an integer as above).
The counts will be stored in it and it will be returned.  The length of
the vector becomes the number of bins.

If the conditions for an external procedure call are satisfied, then it
will be most efficient to make the vector an *INTVEC.

-- Offsetting results in the vector -----------------------------------

It may be useful to place the counts in part of the vector, not
necessarily starting at the first element.  This can be done by passing
a list as __nbins, with three elements:

    _startindex: the index of the first bin
    __nbins: the number of bins
    veclen: the length of the vector

The _hist result will then be of length veclen with the counts in the
elements from _startindex to _startindex + __nbins - 1.

If a vector is to be re-used, it can be given as the third element of
the list, in place of veclen.
*/




