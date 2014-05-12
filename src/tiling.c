#include <R.h>
#include <stdio.h>
#include <Rmath.h>

/* Code for calculating Tiling coefficient.
 * All spike trains must be ordered, smallest time first.
 */

double run_P(int N1, int N2, double dt,
	     double *spike_times_1,
	     double *spike_times_2){
  /* Calculate the term P_1. the fraction of spikes from train 1 that
   * are within +/- dt of train 2.
   */
  
  int i, j, Nab;
		
  Nab=0;
  j=0;
  for(i=0;i<=(N1-1);i++){
    while(j<N2){	
      //check every spike in train 1 to see if there's a spike in
      // train 2 within dt  (don't count spike pairs)
      // don't need to search all j each iteration
      if(fabs(spike_times_1[i]-spike_times_2[j])<=dt){
	Nab=Nab+1;	
	break;				
      }
      else if(spike_times_2[j]>spike_times_1[i]){			
	break;
      }
      else{
	j=j+1;
      }		
    }
  }
  return Nab;
}



double run_T(int N1v, double dtv, double startv, double endv,
	     double *spike_times_1){

  /* Calculate T_A, the fraction of time 'tiled' by spikes with +/- dt.
   *
   * This calculation requires checks to see that (a) you don't count
   * time more than once (when two or more tiles overlap) and checking
   * beg/end of recording.
   */
  
  double dt= dtv;
  double start=startv;
  double end=endv;
  int N1= N1v;
  double time_A;
  int i=0;
  double diff;
  
  //maximum 
  time_A=2*(double)N1*dt;

  // Assume at least one spike in train!
  
  // if just one spike in train 
  if(N1==1){
    
    if((spike_times_1[0]-start)<dt){
      time_A=time_A-start+spike_times_1[0]-dt;
    }
    else if((spike_times_1[0]+dt)>end){
      time_A=time_A-spike_times_1[0]-dt+end;
    }
    
  }
  
  else{				/* more than one spike in train */
    while(i<(N1-1)){
      diff=spike_times_1[i+1]-spike_times_1[i];
      if(diff<2*dt){
	//subtract overlap 	
	time_A=time_A-2*dt+diff;
	
      }
      
      i++;
    }
    
    // check if spikes are within dt of the start and/or end, if so
    // just need to subract overlap of first and/or last spike as all
    // within-train overlaps have been accounted for (in the case that
    // more than one spike is within dt of the start/end
    
    if((spike_times_1[0]-start)<dt){
      time_A=time_A-start+spike_times_1[0]-dt;
    }
    if((end-spike_times_1[N1-1])<dt){
      time_A=time_A-spike_times_1[N1-1]-dt+end;
    }
  }
  return time_A;	
}
	


void run_TM(int *N1v, int *N2v, double *dtv, double *Time,
	    double *index,
	    double *spike_times_1, double *spike_times_2) {

  double TA,  TB,  PA, PB, T;

  int N1= *N1v;
  int N2= *N2v;
  double dt= *dtv;
	
  if(N1==0 || N2==0){
    *index=R_NaN;
  }
  else{
    T=Time[1]-Time[0];
    TA=run_T(N1,dt,Time[0],Time[1], spike_times_1);
    TA=TA/T;
    TB=run_T(N2,dt,Time[0],Time[1], spike_times_2);
    TB=TB/T;
    PA=run_P(N1,N2,dt, spike_times_1, spike_times_2);
    PA=PA/(double)N1;
    PB=run_P(N2,N1,dt, spike_times_2, spike_times_1);
    PB=PB/(double)N2;
    
    *index=0.5*(PA-TB)/(1-TB*PA)+0.5*(PB-TA)/(1-TA*PB);
  }
  
}




void tiling_arr(Sfloat *spikes,
		int *pn,
		int *nspikes,
		int *first_spike,
		Sfloat *rec_time, /* recording time */
		Sfloat *pdt,
		Sfloat *corrs /* return array */) {

  /* Compute all pairwise interactions, include self. */
  /* Elements on lower diagonal are not touched, so those should remain NA. */
  int a, b, n, count;
  Sfloat *sa, *sb; 		/* pointers to current spike trains  */
  int n1, n2;
  Sfloat k1, res;
  int debug;

  n = *pn;
	    
  for (a=0; a<n; a++) {
    n1 = nspikes[a];
    sa = &(spikes[first_spike[a]]);
    
    for (b=a; b<n; b++) {
      n2 = nspikes[b];
      sb = &(spikes[first_spike[b]]);

      run_TM(&n1, &n2, pdt, rec_time, &res, sa, sb);
      corrs[(b*n)+a] = res;
      
    }
  }
}
