/***************************************************
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/

#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <time.h>

// #######  Fundamental constants  #######
#define hbar ( 0.658212196 )     // eVfs // meVps // �eVns
#define PI M_PI

// #######  Definition of time domain  #######
#define MAX_STEPS (8000)
double delta_t = 0.5 ;


// #######  Detector resolution  #######
#define RESOLUTION ( 100.000 / (2.355) ) // ps



// #######  Spectral properties  #######
//#define FSS         (    +0.0115 / hbar  ) // in meV --> ps^-1
// Parameter:
// Gershoni:  FSS = 0.052 * hbar = 34,0 ?eV
// Samir QD1: FSS = 0.005 * hbar =  3,5 ?eV
// Samir QD2: FSS = 0.017 * hbar = 11,5 ?eV
// Samir QD3: FSS = 0.071 * hbar = 47,0 ?eV

#define DELTA_XX    ( -1.000 / hbar ) // Biexciton shift

// #######  Lindblad Parameters  #######

// Radiative Recombination
#define GAM_R    ( 0.0024 * 1.0 )
#define GAM_R_XH ( GAM_R * 1.00 ) // 1/ps --> 1/T1 = 500ps
#define GAM_R_XV ( GAM_R * 1.00 )
#define GAM_R_BH ( GAM_R * 1.00 )
#define GAM_R_BV ( GAM_R * 1.00 )

// Pure Dephasing
#define GAM_P ( 0.017 )
#define GAM_PH ( 0.017 )
#define GAM_PV ( 0.017 )

#define GAM_P_XH ( GAM_PH  * 1. * 1.0 ) // 1/ps -- 1/T2* = 50ps
#define GAM_P_XV ( GAM_PV  * 1. * 1.0 ) // 1/ps -- 1/T2* = 50ps
#define GAM_P_BH ( GAM_PH  * 1. * 1.0 ) // 1/ps -- 1/T2* = 50ps
#define GAM_P_BV ( GAM_PV  * 1. * 1.0 ) // 1/ps -- 1/T2* = 50ps
//------------------------------------------------------------




// #########################################################
// ############# RK4 #######################################
// #########################################################

#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})

#define SECHSTEL 0.1666666666666666666666667
#define DRITTEL  0.3333333333333333333333333
#define HALFTE   0.5
