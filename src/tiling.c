#include <R.h>
#include <time.h>
#include <stdio.h>
#include <Rmath.h>

double run_P(int N1,int N2,double dt,double *spike_times_1,double *spike_times_2){

	int i;
	int j;
	int Nab;
		
	Nab=0;
	j=0;
	for(i=0;i<=(N1-1);i++){
		while(j<N2){	
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



double run_T(int N1v,double dtv,double Tv,double *spike_times_1,double *spike_times_2){
	
	int i;
	int j;
	int u;
	double dt= dtv;
	double start;
	double end;
	int N1= N1v;
	double T= Tv;
	double time_A;
	
	if((spike_times_1[1]-spike_times_1[0])> 2*dt){
		if(spike_times_1[0]>dt){
			time_A=2*dt;				
		}
		else{					
			time_A=spike_times_1[0]+dt;
		}
		u=1;
	}
	else{
		start=spike_times_1[0];
		end=spike_times_1[1];
		for(j=0;j<(N1-1);j++){
			if((spike_times_1[j+1]-spike_times_1[j])<=2*dt){
				end=spike_times_1[j+1];
			}
			else{
				break;
			}		
		}
		if(spike_times_1[0]<=dt){
			time_A=end+dt;
		}
		else{
			time_A=end-start+2*dt;
		}
		u=j+1;
	}

	i=u;

//n.b. if second spike is within dt of first spike, it'll be picked up by the first loop
	
	while(i <=(N1-1)){
		if(i==(N1-1)){
			if((T-spike_times_1[i])<=dt){				
				time_A=time_A+dt+T-spike_times_1[i];
			}
			else{
				time_A=time_A+2*dt;				
			}
			i++;
		}else{
			if((spike_times_1[i+1]-spike_times_1[i])<=2*dt){
				start=spike_times_1[i];
				end=spike_times_1[i+1];
				
				for(j=i;j<(N1-1);j++){
					
					if((spike_times_1[j+1]-spike_times_1[j])<=2*dt){
						end=spike_times_1[j+1];	
					}
						else{
							break;
					}
					
				}
				if((end+dt)>T){
					time_A=time_A+(T-start)+dt;					
				}
				else{
					time_A=time_A+(end-start)+2*dt;						
				}
				i=j+1;	
			}
			else{
				time_A=time_A+2*dt;
				i=i+1;
			}
		}
	}
	return time_A;
}

void run_TM(int *N1v,int *N2v,double *dtv,double *Tv,double *index,double *spike_times_1,double *spike_times_2){

	double TA;
	double TB;
	double PA;
	double PB;
	int N1= *N1v;
	int N2= *N2v;
	double dt= *dtv;
	double T=*Tv;

	TA=run_T(N1,dt,T, spike_times_1, spike_times_2);
	TA=TA/T;
	TB=run_T(N2,dt,T, spike_times_2, spike_times_1);
	TB=TB/T;
	PA=run_P(N1,N2,dt, spike_times_1, spike_times_2);
	PA=PA/(double)N1;
	PB=run_P(N2,N1,dt, spike_times_2, spike_times_1);
	PB=PB/(double)N2;
	
	*index=0.5*(PA-TB)/(1-TB*PA)+0.5*(PB-TA)/(1-TA*PB);

}




