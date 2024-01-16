/*************************************************** 
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/
// C++ Interface: header
// Description: 
// Author: Alexander Carmele <alex@itp.tu-berlin.de>, (C) 2009
// Copyright: See COPYING file that comes with this distribution

#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
 
#define hbar 0.658212196      // eVfs 

// MPS code time time 300/dt = 600, however I subtract a tau intervall
// time input 300 --> 600 steps, feedback time 100 --> 500 steps * dt = t --> 250 
// 300/dt-100 feedback_time , e.g. for dt 0.5 --> 600-100 = 500 

#define MPS_time_end (300)
#define MPS_dt (0.3)
#define MPS_feedback_time  (50-1)
#define MPS_phi ( -10.35  )

double TIME_END=MPS_time_end -MPS_dt*MPS_feedback_time;
double FEEDBACK_TIME = MPS_dt*MPS_feedback_time;
int N=17000;  
//double delta_t = 0.1;// TIME_END/N;

// #########################################################
// ########## PHOTON PARAMETER #############################
// #########################################################
 #define KAPPA   ( 0.1 * 0.1 ) 
 
// #########################################################
// ########################## RK4 ##########################
// #########################################################


#define ce       ( *(derivates     + 1 )) 
#define cg       ( *(derivates     + 2 )) 
 
#define ce_OUT       ( *(derivates_out + 1 )) 
#define cg_OUT       ( *(derivates_out + 2 )) 

 int N_DGL=4;

// ############# RK4 ##########################
#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})
// ############################################


//////////////////////////////////////////////////////////////////////////////
// ########################## DGL-System #####################################
//////////////////////////////////////////////////////////////////////////////


void calculate_next_time_step (complex double* derivates,complex double* derivates_out,double t,complex double* memory_cg,int now,int N_tau)
{   
   
if ( now > N_tau) cg_OUT = - KAPPA * cg + cexp(I*MPS_phi)*KAPPA * memory_cg[now - N_tau]; 
else              cg_OUT = - KAPPA * cg;

return ;
}

// ###################### Ende     ############################################  
///////////////////////////////////////////////////////////////////////////////  
// #########################  MAIN    ######################################### 

int main()
{
    
double delta_t = TIME_END/N;    
int N_tau=(int)(FEEDBACK_TIME/delta_t); 
printf("N_tau=%i \n",N_tau);    
    
int i;
complex double *derivates = calloc(2*N_DGL,sizeof(double));
complex double *temp = calloc(2*N_DGL,sizeof(double));
complex double *k1 = calloc(2*N_DGL,sizeof(double));
complex double *k2 = calloc(2*N_DGL,sizeof(double));
complex double *memory_cg = calloc(2*N,sizeof(double));

cg = 1.0;
  
FILE *f_dat;
f_dat = fopen("RK4_FEEDBACK_SE.dat","w");
double t;
for (int k=0; k<N; k++)
{  
 t=k*delta_t;   
 memory_cg[k] = cg;
 fprintf(f_dat,"%.10f \t %.10f \n",delta_t*k,creal(cg*conj(cg)));

 calculate_next_time_step(derivates, k1 , t, memory_cg, k,N_tau); 
 CALCULATE_TEMP(derivates,k1,delta_t*0.5,temp,i);
 calculate_next_time_step(temp, k2, t+delta_t*0.5, memory_cg, k,N_tau); 
 CALCULATE_TEMP(derivates,k2,delta_t*0.5,temp,i);
 ADD_VECTOR(k1,k2,2.0,i);
 calculate_next_time_step(temp, k2, t+delta_t*0.5, memory_cg, k,N_tau); 
 CALCULATE_TEMP(derivates,k2,delta_t,temp,i);
 ADD_VECTOR(k1,k2,2.0,i);
 calculate_next_time_step(temp, k2, t+delta_t, memory_cg, k,N_tau);
 ADD_VECTOR(k1,k2,1.0,i);
 ADD_VECTOR(derivates,k1,delta_t*0.166666666666667,i);
}
    
free(temp); free(derivates); free(k1); free(k2);
free(memory_cg);
fclose(f_dat);
return 0;
}

