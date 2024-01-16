/***************************************************
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/

 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>
 #include <time.h>

// #define DOUBLE_TYPE double

  #define hbar 0.658212196      // eVfs
 #define PI 3.14159265358979323846264338327

 #define TIME_DYNAMICS_TO_NS 0.0175600
 #define T_STEPS 5400
 int STEP_IS_SAVED= 1; // jeder ... Schritt wird in ein File geschrieben!!
 int N=T_STEPS;
 double delta_t = ((double)TIME_DYNAMICS_TO_NS) * 1000. * 1000. /((double)T_STEPS);


 #define c_0 299.8  // nm/fs
 // ##### Beachte, die Zeit skaliert mit 2 PI wegen exp(I*t*2*PI/T) #######
 #define FEEDBACK_TIME_TAU ( 0.025 * 0.1 * PI) //( 0.075 * 2.33 ) //ns
 #define L                 ( FEEDBACK_TIME_TAU * 1000.0 * c_0 * 0.5 )
 #define omega_0 1.500   // 1/fs
 #define G_fak (0.15*1.25*1.) // feedback strength
 #define WITH_MIRROR
 // k-bereich muss je laenger die zeitdynamik ist, um so breiter gewaehlt werden
 #define DISCRETIZE_K_TO 0.0051
 #define K_INITIAL 0.0049
 #define N_k 500
 double delta_k = (DISCRETIZE_K_TO-K_INITIAL)/N_k;
 double M = (0.00008*0.); // (fs)^(-1);
// bspw. g_fak 0.2, with_mirror 0.00498 - 0.00502, 0.15ns, t_steps 305, N_k 300
// ein delta_k von mindestens 1e-07 wichtig
// ... k * c_0, also 0.00495 --> 1.48401, Intervall: 0.00505 --> 1.51399 1/fs


// #########################################################
// ########################## RK4 ##########################
// #########################################################

#define ce        ( *(derivates + 1 ))
#define cg        ( *(derivates + 2 ))
#define ck(p)     ( *(derivates + 3 + (p) ))

#define ce_OUT        ( *(derivates_out + 1 ))
#define cg_OUT        ( *(derivates_out + 2 ))
#define ck_OUT(p)     ( *(derivates_out + 3 + (p) ))

#define Zahl_der_DGL ( 15 + 1 * (N_k + 1)  )
int N_DGL=Zahl_der_DGL;

// ############# RK4 ##########################
#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})




//////////////////////////////////////////////////////////////////////////////
// ########################## DGL-System #####################################
//////////////////////////////////////////////////////////////////////////////


void calculate_next_time_step
(
complex double *derivates,
complex double *derivates_out,
        double t,
complex double *rot_fak,
        double *Gk
)
{

int p;

ce_OUT   = + I * M * cg; // cavity dynamics
cg_OUT   = + I * M * ce;

for (p=0;p<N_k;p++)
{
  cg_OUT   += - I * Gk[p]  * ck(p) * rot_fak[p]  * delta_k;
  ck_OUT(p) = - I * conj(Gk[p])  * cg  * conj(rot_fak[p]);
}
return ;
}


 int main()
 {  int k,i;
    double t;
    complex double *derivates = calloc(2*N_DGL,sizeof(double));
   // Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
    complex double *temp = calloc(2*N_DGL,sizeof(double));
    // Die errechneten Zwischenfunktionswerte - die Steigungen
    complex double *k1 = calloc(2*N_DGL,sizeof(double));
    complex double *k2 = calloc(2*N_DGL,sizeof(double));


FILE *f_dat;
f_dat=fopen("feedback_modes_original.dat","w");
//////////////////////////////////////////////////////////////////////////////////
// ###############################################################################
// ############################## INITIAL CONDITIONS #############################
// ###############################################################################

  ce = 0.0;
  cg = 1.0;

// ###########################################################################
// ########################## TIME _DOMAIN SOLUTION ##########################
// ###########################################################################

#ifdef WITH_MIRROR
printf("=============================\n");
printf("Delay tau=%.10f ns\n",0.002*L/c_0);
printf("Distance L=%.1f nm ----> L=%.1f mm\n",L*1000.,L/1000.);
printf("=============================\n");
#else
printf("====== No Delay ========\n");
#endif

// ####################################################
// ########## optimize : rf factor - Phononen #########
// ####################################################

double temp_om = 0.;
complex double *rf_fak_1    = calloc(2*N_k,sizeof(double));
complex double *rf_fak_2    = calloc(2*N_k,sizeof(double));
complex double *rf_fak_dt   = calloc(2*N_k,sizeof(double));

for (i=0; i<N_k; i++)
 {
  temp_om = omega_0 - (K_INITIAL + delta_k * i) * c_0;
  rf_fak_dt[i] = cexp( - I * temp_om * delta_t * 0.5 );
  rf_fak_1[i]  = 1. + I * 0.;
  rf_fak_2[i]  = 1. + I * 0.;
 }
// #########################################
// ########## END: optimize : rf factor ####
// #########################################

// #################################################
// ###### Begin: Optimize Mirror coupling ##########
// #################################################

double *coupl_fak = calloc(N_k,sizeof(double));

for (i=0; i<N_k; i++)
{
  #ifdef WITH_MIRROR
  coupl_fak[i] = G_fak * sin( ( K_INITIAL + delta_k * i ) * L * 1000. );
  #else
  coupl_fak[i] = G_fak;
  #endif
}
// #################################################
// ###### End: Optimize Mirror coupling ##########
// #################################################
//Beginn der Integration von 0fs bis N*Schrittweite in fs
   for (k=0; k<N; k++)
     {
       t=delta_t*k;

      if ( k % STEP_IS_SAVED ==0) { fprintf(f_dat,"%.10f \t %.10f \t %.10f \n",delta_t*k*0.000001,creal(cg*conj(cg)),creal(ce*conj(ce))); }
// ################################################################
// ############### 4 VECTOR RUNGE KUTTA ###########################
// ################################################################
          calculate_next_time_step(derivates, k1 , t, rf_fak_2 , coupl_fak);
          CALCULATE_TEMP(derivates,k1,delta_t*0.5,temp,i);

          // #### RF Photon ####################################
          for (i=0; i<N_k; i++)
          {
           rf_fak_1[i] = rf_fak_2[i] * rf_fak_dt[i]; //berechne
           rf_fak_2[i] = rf_fak_1[i] * rf_fak_dt[i];
    	  }

      calculate_next_time_step(temp, k2, t+delta_t*0.5, rf_fak_1, coupl_fak);
	  CALCULATE_TEMP(derivates,k2,delta_t*0.5,temp,i);
	  ADD_VECTOR(k1,k2,2.0,i);

      calculate_next_time_step(temp, k2, t+delta_t*0.5, rf_fak_1, coupl_fak);
	  CALCULATE_TEMP(derivates,k2,delta_t,temp,i);
  	  ADD_VECTOR(k1,k2,2.0,i);

      calculate_next_time_step(temp, k2, t+delta_t, rf_fak_2, coupl_fak);
	  ADD_VECTOR(k1,k2,1.0,i);
	  ADD_VECTOR(derivates,k1,delta_t*0.166666667,i);
// ################################################################
// ############ END OF 4 VECTOR RUNGE KUTTA #######################
// ################################################################

    }

fclose(f_dat);

  free(temp);
  free(derivates);
  free(k1);
  free(k2);
  free(rf_fak_1);
  free(rf_fak_2);

  free(coupl_fak);

  return 0;
}
