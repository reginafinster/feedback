/*************************************************** 
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/
// C++ Interface: header
// Description: 
// Author: Alexander Carmele <alex@itp.tu-berlin.de>, (C) 2009
// Copyright: See COPYING file that comes with this distribution

#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <time.h>
 
#define DOUBLE_TYPE double

#define hbar 0.658212196      // eVfs 
#define PI M_PI

// ################ ZEITRAUM #######################

#define T_STEPS 600

    int N=T_STEPS;  

 DOUBLE_TYPE delta_t = 0.1;

// #########################################################
// ########## PHOTON PARAMETER #############################
// #########################################################

 #define c_0 299.8  // nm/fs
 // ##### Beachte, die Zeit skaliert mit 2 PI wegen exp(I*t*2*PI/T) #######
 #define FEEDBACK_TIME_TAU ( 0.025 * 1. * PI) //( 0.075 * 2.33 ) //ns
 #define N_tau   ( 98 )
 #define L       ( FEEDBACK_TIME_TAU * 1000.0 * c_0 * 0.5 )  
 #define KAPPA   ( 1.0 ) 
 
// #########################################################
// ########################## RK4 ##########################
// #########################################################


#define ce       ( *(derivates     + 1 )) 
#define cg       ( *(derivates     + 2 )) 
 
#define ce_OUT       ( *(derivates_out + 1 )) 
#define cg_OUT       ( *(derivates_out + 2 )) 

 
#define Zahl_der_DGL ( 15 + 100 )

int N_DGL=Zahl_der_DGL;

// ############# RK4 ##########################
#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})

#define SECHSTEL 0.1666666666666666666666667
#define DRITTEL  0.3333333333333333333333333
#define HALFTE   0.5
// ############################################


//////////////////////////////////////////////////////////////////////////////
// ########################## DGL-System #####################################
//////////////////////////////////////////////////////////////////////////////


void calculate_next_time_step 
(
complex DOUBLE_TYPE* derivates, 
complex DOUBLE_TYPE* derivates_out,  
        DOUBLE_TYPE t, 
complex DOUBLE_TYPE* memory_cg,
        int N
)
{   

complex DOUBLE_TYPE temp;
 temp = 0.;

// TODO Annahme: im MPS verdoppelt sich Gamma ab t>tau	   
if ( N > N_tau) cg_OUT = - 2.28*KAPPA * cg + 2.28*KAPPA * memory_cg[N - N_tau]; 
else 
    cg_OUT   = - 1.0*KAPPA * cg;

return ;
    }

// ###################### Ende     ############################################  
///////////////////////////////////////////////////////////////////////////////  
// #########################  MAIN    ######################################### 

 int main()                                     
 {  int k,n,o;
    DOUBLE_TYPE q,t;
    double zeit_start,zeit_start_2,zeit_ende, v,zwischen_zeit,rest_zeit;
    DOUBLE_TYPE complex temp_z;
    DOUBLE_TYPE temp1,temp2,temp3,temp4,temp5,temp6,temp7, temp0,temp8,temp9,temp10,temp11,temp12;
    complex DOUBLE_TYPE *derivates = calloc(2*N_DGL,sizeof(DOUBLE_TYPE));
   // Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
    complex DOUBLE_TYPE *temp = calloc(2*N_DGL,sizeof(DOUBLE_TYPE));
    // Die errechneten Zwischenfunktionswerte - die Steigungen
    complex DOUBLE_TYPE *k1 = calloc(2*N_DGL,sizeof(DOUBLE_TYPE));
    complex DOUBLE_TYPE *k2 = calloc(2*N_DGL,sizeof(DOUBLE_TYPE));

   temp1 = 0.; 
   temp2 = 0.;

// ###############################################################################
// ############################## INITIAL CONDITIONS #############################
// ###############################################################################

int i,l;

  cg = 1.0;
  
// ###########################################################################
// ########################## TIME _DOMAIN SOLUTION ##########################
// ###########################################################################

FILE *f_dat;
f_dat = fopen("RK4_FEEDBACK_SE.dat","w");

// ####################################################
// ########## optimize : rf factor - Phononen #########
// ####################################################


complex DOUBLE_TYPE *memory_cg = calloc(2*N,sizeof(DOUBLE_TYPE));

DOUBLE_TYPE DT=delta_t;
//Beginn der Integration von 0fs bis N*Schrittweite in fs                            
   for (k=0; k<N; k++)
     {  
       t=delta_t*k;
       // #################################
       temp_z=cg;
       temp2 =cg*conj(cg);
       memory_cg[k] = temp_z;
       // ###############################
       fprintf(f_dat,"%.10f \t %.10f \n",t,temp2);
       if (temp1 <5.00) 
      {
// ################################################################
// ############### 4 VECTOR RUNGE KUTTA ###########################
// ################################################################
      calculate_next_time_step(derivates, k1 , t, memory_cg, k); 
      CALCULATE_TEMP(derivates,k1,DT*HALFTE,temp,i);
      calculate_next_time_step(temp, k2, t+delta_t*0.5, memory_cg, k); 
	  CALCULATE_TEMP(derivates,k2,DT*HALFTE,temp,i);
	  ADD_VECTOR(k1,k2,2.0,i);
      calculate_next_time_step(temp, k2, t+delta_t*0.5, memory_cg, k); 
	  CALCULATE_TEMP(derivates,k2,DT,temp,i);
  	  ADD_VECTOR(k1,k2,2.0,i);
      calculate_next_time_step(temp, k2, t+delta_t, memory_cg, k);
	  ADD_VECTOR(k1,k2,1.0,i);
	  ADD_VECTOR(derivates,k1,DT*SECHSTEL,i);
// ################################################################
// ############ END OF 4 VECTOR RUNGE KUTTA #######################
// ################################################################
        }
    }
    
  free(temp);
  free(derivates);
  free(k1);
  free(k2);
  free(memory_cg);
  fclose(f_dat);

  return 0;
}

