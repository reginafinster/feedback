//Programm kann entweder steady state suchen oder eine vorgegebene laufzeit suchen, je nachdem was strenger ist (MAX_DEVIA oder TIME_DYN...). 

/*************************************************** 
 *   Copyright (C) 2008 by Alexander Carmele   *
 *   221b Baker St, Marylebone, London W1.     *
 ***************************************************/

 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>
 #include <time.h>
 
 #define DOUBLE_TYPE double
// C++ Interface: header
// Description: 
// Author: Alexander Carmele <alex@itp.tu-berlin.de>, (C) 2009
// Copyright: See COPYING file that comes with this distribution

#define hbar 0.658212196      // eVfs 
#define PI M_PI

// ################ ZEITRAUM #######################
#define TIME_DYNAMICS_TO_NS 0.00100000
#define T_STEPS 5000
int N = T_STEPS;
int STEP_IS_SAVED= 1; // jeder ... Schritt wird in ein File geschrieben!!
int TIME_OUTPUT  = 12; //nur für terminal print
#define N_OUT ( 1000 ) //groesse des ergebnis arrays - vorsicht, stst auf 999

#define MAX_DEVIA 0.0000000001 //kriterium für steady state

#define J_CONST ( 0.1  * 1.0    ) // value for spin chain expressed in S_x/y --> 0.5 * sigma_x/y
#define V_CONST ( 0.1  * 1.0    ) // spin half picture!!!
#define H_CONST ( 0. * 1.00 )  // spin half picture!!

//#define OWN_DISORDER	//erzeugt selbst disorder für abgleich der gemittelten werte
#define NO_AVERAGE 10
//#define H_MAX 1.0 //entweder hier festlegen oder auf argv[1] übergeben (castet dann ein double aus dem char, hoffentlich :)


//#define SHIFT_IMPORT //importiert shifts aus mps für abgleich einzelner runs, ausmachen für OWN_DISORDER
#define DISFACTOR 2. //mit -2. passt es zu Reginas MPS, Werte werden aus shifts.dat eingelesen
	
#define N_spins ( 4 )

//#define INITIAL_STEADY
//#define INITIAL_NEEL_0_UP //Neel-state als anfangszustand anstatt GZS, 1st site up
//#define INITIAL_NEEL_0_DOWN //1st site down
//#define INITIAL_ALLUP 
#define ALL_DOWN_LAST_UP

DOUBLE_TYPE delta_t = ((DOUBLE_TYPE)TIME_DYNAMICS_TO_NS) * 1000. * 1000. /((DOUBLE_TYPE)T_STEPS); 

#define RandPrec (10000)

#define V(i)     ( *(delta_array       + (i) ) )
#define J(i)     ( *(coupling_array    + (i) ) )
#define H(i)     ( *(magnetic_field    + (i) ) )

#define MUE		 (0.00) //für REPRO LEON 
#define GAMMA	 (0.0)

#define LOSS     (0.24*0.24) // für BM FB-Programm: Gamma von dort quadrieren!!
//#define LOSS     (0.5*2.2*2.2) // sigma_plus |0><1| jump operator
//#define LOSS     (GAMMA*(1.+MUE)) // sigma_plus |0><1| jump operator
#define PUMP     (GAMMA*(1.+MUE)) // sigma_minus |1><0| jump operator

#define LOSS_REV (GAMMA*(1.-MUE)) //für REPRO LEON
#define PUMP_REV (GAMMA*(1.-MUE))

#define MIDDLE_POS (4)
#define LOSS_MIDDLE (0.0)
#define PUMP_MIDDLE (0.0)

//#define LONGRANGE	//schaltet Longrange Term ein.
//#define alpha 1000.	//longrange parameter bei Leon
//double A;			//Normierung (wird in function set_A befüllt)S

// #########################################################
// ############# RK4 #######################################
// #########################################################

#define CALCULATE_TEMP(a,b,cb,dest,i) ({for(i=0;i<N_DGL;i++) (dest)[i]=(a)[i]+(cb)*(b)[i];})
#define ADD_VECTOR(dest,in,in_fak,i) ({for(i=0;i<N_DGL;i++) dest[i]=dest[i]+(in_fak)* (in)[i];})

#define SECHSTEL 0.1666666666666666666666667
#define DRITTEL  0.3333333333333333333333333
#define HALFTE   0.5

#define NORM ( 1./(2.*(N_spins-1.)) ) //für current-berechnung

#define PLOT //plottet mit gnuplot
#define CURRENTPLOT //schreibt einzelne <j>s mit

// #################################################################################
// ############## DGL system  ######################################################
// #################################################################################

void STEADY_STATE_TEST ( 
complex double *derivates,
complex double *derivates_stst,
		double *OUTPUT,
		int N_dgl
)
{
 int x;
 OUTPUT[999]=0.;
 for (x=0;x<N_dgl;x++)
 {
   OUTPUT[999] += cabs(derivates[x]-derivates_stst[x]);
   derivates_stst[x] = derivates[x];
   }
 //printf("StstTest = %f \n",OUTPUT[999]);
 return;
}

#ifdef LONGRANGE
double Set_A(){
	double A=0.0;
	for(int i=0;i<N_spins-1;i++){ for(int j=i+1;j<N;j++) A += 1/((N_spins-1)*pow(cabs(j-i),alpha));}
	printf("A %f\n",A);
	return A;
}
#endif

void calculate_next_time_step 
(
complex double *derivates, 
complex double *derivates_out, 
        double t,
        double *delta_array, 
        double *coupling_array,
        double *magnetic_field,
        int    *POW_2,
        double *OUTPUT
)
{   

int k;
int i;

int POS_TARGET;
int POS_SOURCE;
int L[N_spins+1];
int R[N_spins+1];
int TRACE_SWITCH;
int CURRENT_SWITCH;
int FINDER;
int REST;
double s_z_temp=0.;

for(i=0;i<N_OUT;i++) OUTPUT[i]=0.; // initialisiere ergebnisarray

// ----------- linke seite dichtematrix ----------- 

for (k=0;k<POW_2[2*N_spins];k++) 
{ 
POS_TARGET=k;
TRACE_SWITCH = 0; // diagonal element or not
CURRENT_SWITCH = 0; // how far off diagonal? 2 -> potential candidate for the current
// ---------- attribute which site is occupied which not depending on state 'k' ---------------------------
// we start from left <123 ... N_spins|rho|123 ... N_spins>, 
FINDER = POS_TARGET;

for (i=0;i<N_spins;i++) 
  {
  REST = FINDER % POW_2[i+1];    
  L[i] = REST/POW_2[i]; 
  FINDER = FINDER - REST;
  }
for (i=0;i<N_spins;i++) 
  {
  REST = FINDER % POW_2[N_spins+i+1];    
  R[i] = REST/POW_2[N_spins+i]; 
  FINDER = FINDER - REST;
  if ( L[i] != R[i] )  {TRACE_SWITCH = -1; CURRENT_SWITCH += 1;}
  }
// now we have for every state betweek 0 - 2^(2_N_states-1) whether the site is occupied or not
// k=POS_TARGET=L[0]*POW_2[0]   +L[1]*POW_2[1]     +...+L[N_spins-1]*POW_2[N_spins-1]
//             +R[0]*POW_2[N_spins]+R[2]*POW_2[N_spins+1]+...+R[N_spins-1]*POW_2[2*N_spins-1]

//printf("k=%i -- ",k); for(i=0;i<N_spins;i++) printf("L[%i]=%i --",i,L[i]); for(i=0;i<N_spins;i++) printf("R[%i]=%i --",i,R[i]); printf("\n");

derivates_out[POS_TARGET] = 0.; //initialize 
for (i=0;i<N_spins;i++)
{
// conventionally sigma_z = |0><0| - |1><1|, sigma_plus = |0><1|    
// -------------------------------------------------------------------------------------------------------
// ----------------------------------------- LEFT - SIDE - DENSITY MATRIX (globally minus sign) ----------
// -------------------------------------------------------------------------------------------------------

// single site Hamiltonian contribution    
if (L[i]==1)  { derivates_out[POS_TARGET] += -I * (-1.*0.5)* H(i)  *derivates[POS_TARGET]; }
if (L[i]==0)  { derivates_out[POS_TARGET] += -I * ( 1.*0.5)* H(i)  *derivates[POS_TARGET];}
// end of single site Hamiltonian contribution

// interaction Hamiltonian contribution (nearest neighbor)
// remember V interaction if the sites are in the same state Hamiltonian sign plus, if not minus
// remember J interaction only if sites not in the same state Hamiltonian sign plus
if (i<N_spins-1) // last site does not interact with first site, open boundary
{ 
if (L[i]==1) 
   {
   if ( L[i+1]==1 )  
      {
       derivates_out[POS_TARGET]+= -I * ( 1.) * V(i)     * derivates[POS_TARGET] ;
      }
   if ( (L[i+1]==0) ) 
      {//L[i]==1,L[i+1]==0 is driven via states with L[i]==0, and L[i+1]==1
       POS_SOURCE = POS_TARGET + (0-1)*POW_2[i] + (1-0)*POW_2[i+1] ;   
       derivates_out[POS_TARGET] += -I * ( 2.) * J(i)     * derivates[POS_SOURCE]   
                                    -I * (-1.) * V(i)     * derivates[POS_TARGET] ;
       } 
   }
if (L[i]==0)  
   {
   if (L[i+1]==1)  
      { //L[i]==0,L[i+1]==1 is driven via states with L[i]==1, and L[i+1]==0
       POS_SOURCE = POS_TARGET + POW_2[i] - POW_2[i+1]; 
       derivates_out[POS_TARGET]+= -I * ( 2.) * J(i)     * derivates[POS_SOURCE]  
                                   -I * (-1.) * V(i)     * derivates[POS_TARGET] ;
      }
   if (L[i+1]==0)  
      {       
       derivates_out[POS_TARGET]+= -I * ( 1.) * V(i)     * derivates[POS_TARGET] ;
      }
    }
} 
// end of interaction Hamiltonian contribution (nearest neighbor)
// --------------------------------------------------------------------------------------------------------
// ----------------------------------------- RIGHT - SIDE - DENSITY MATRIX --------------------------------
// --------------------------------------------------------------------------------------------------------

// single site Hamiltonian contribution    
if (R[i]==1)  { derivates_out[POS_TARGET] += +I * (-1.*0.5) * H(i)  *derivates[POS_TARGET]; }
if (R[i]==0)  { derivates_out[POS_TARGET] += +I * ( 1.*0.5) * H(i)  *derivates[POS_TARGET]; }
// end of single site Hamiltonian contribution

// interaction Hamiltonian contribution (nearest neighbor)
if (i<N_spins-1) // last site does not interact with first site, open boundary
{
if (R[i]==1)   
   {
    if (R[i+1]==1)  
       { 
        derivates_out[POS_TARGET]+= +I * ( 1.) * V(i) * derivates[POS_TARGET] ; 
       } 
    if (R[i+1]==0) 
       {
        POS_SOURCE = POS_TARGET - POW_2[N_spins+i] + POW_2[N_spins+i+1]  ; 
        derivates_out[POS_TARGET]+= +I * ( 2.) * J(i)     * derivates[POS_SOURCE]  
                                    +I * (-1.) * V(i)     * derivates[POS_TARGET] ;
       }
   }
if (R[i]==0)   
   { 
    if ( R[i+1]==1 ) 
       { 
        POS_SOURCE = POS_TARGET + POW_2[N_spins+i] - POW_2[N_spins+i+1]; 
        derivates_out[POS_TARGET]+= +I * ( 2.) * J(i)     * derivates[POS_SOURCE]  
                                    +I * (-1.) * V(i)     * derivates[POS_TARGET] ;
       }
    if (R[i+1]==0) 
       { 
        derivates_out[POS_TARGET]+= +I * ( 1.) * V(i)     * derivates[POS_TARGET] ;
       }
   }
} 


#ifdef LONGRANGE
//Long-Range (LEON)
for(int j=i+1;j<N_spins;j++){
if( (L[i]!=L[j]) && (R[i]==R[j]) ){		derivates_out[POS_TARGET] += I*0.5*J_CONST/(A*pow(j-i,alpha))*derivates[POS_TARGET];}
if( (L[i]==L[j]) && (R[i]!=R[j]) ){		derivates_out[POS_TARGET] -= I*0.5*J_CONST/(A*pow(j-i,alpha))*derivates[POS_TARGET];}
}
#endif
  
// --------------------------------------------------------------------------------------------------
// ------------------------------------- LINDBLAD ---------------------------------------------------
// --------------------------------------------------------------------------------------------------
 if(i==0){//ohne diese schleifen N_spins mal pro durchlauf!! daher kam der faktor.
 if   (L[0]==0)                 derivates_out[POS_TARGET] += -1.*PUMP*derivates[POS_TARGET];
 if   (R[0]==0)                 derivates_out[POS_TARGET] += -1.*PUMP*derivates[POS_TARGET];
 if ( (L[0]==1) && (R[0]==1) ) {
                                POS_SOURCE = POS_TARGET - POW_2[0] - POW_2[N_spins]; 
                                derivates_out[POS_TARGET] += +2.*PUMP*derivates[POS_SOURCE];
                               }
 }
 if(i==N_spins-1){
 if   (L[N_spins-1]==1)            derivates_out[POS_TARGET] += -1.*LOSS*derivates[POS_TARGET];
 if   (R[N_spins-1]==1)            derivates_out[POS_TARGET] += -1.*LOSS*derivates[POS_TARGET];
 if ( (L[N_spins-1]==0) && (R[N_spins-1]==0)) 
                               {
                                POS_SOURCE = POS_TARGET + POW_2[N_spins-1] + POW_2[N_spins+N_spins-1]; 
                                derivates_out[POS_TARGET] += +2.*LOSS*derivates[POS_SOURCE];
                               }
 }
// --------------------------------------------------------------------------------------------------
// ------------------------------------- EOF LINDBLAD ---------------------------------------------------
// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
// ------------------------------------- LINDBLAD REVERSE FOR LEON REPRO---------------------------------
// --------------------------------------------------------------------------------------------------
 if(i==N_spins-1){
 if   (L[N_spins-1]==0)                 derivates_out[POS_TARGET] += -1.*PUMP_REV*derivates[POS_TARGET];
 if   (R[N_spins-1]==0)                 derivates_out[POS_TARGET] += -1.*PUMP_REV*derivates[POS_TARGET];
 if ( (L[N_spins-1]==1) && (R[N_spins-1]==1) ) {
                                POS_SOURCE = POS_TARGET - POW_2[N_spins-1] - POW_2[N_spins+N_spins-1]; 
                                derivates_out[POS_TARGET] += +2.*PUMP_REV*derivates[POS_SOURCE];
                               }
 }
 if(i==0){
 if   (L[0]==1)            derivates_out[POS_TARGET] += -1.*LOSS_REV*derivates[POS_TARGET];
 if   (R[0]==1)            derivates_out[POS_TARGET] += -1.*LOSS_REV*derivates[POS_TARGET];
 if ( (L[0]==0) && (R[0]==0)) 
                               {
                                POS_SOURCE = POS_TARGET + POW_2[0] + POW_2[N_spins]; 
                                derivates_out[POS_TARGET] += +2.*LOSS_REV*derivates[POS_SOURCE];
                               }
 }
// --------------------------------------------------------------------------------------------------
// ------------------------------------- EOF LINDBLAD REVERSE---------------------------------------------------
// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
// ------------------------------------- LINDBLAD MIDDLE_POS---------------------------------
// --------------------------------------------------------------------------------------------------
 if(i==MIDDLE_POS){
 if   (L[MIDDLE_POS]==0)                 derivates_out[POS_TARGET] += -1.*PUMP_MIDDLE*derivates[POS_TARGET];
 if   (R[MIDDLE_POS]==0)                 derivates_out[POS_TARGET] += -1.*PUMP_MIDDLE*derivates[POS_TARGET];
 if ( (L[MIDDLE_POS]==1) && (R[MIDDLE_POS]==1) ) {
                                POS_SOURCE = POS_TARGET - POW_2[MIDDLE_POS] - POW_2[N_spins+MIDDLE_POS]; 
                                derivates_out[POS_TARGET] += +2.*PUMP_MIDDLE*derivates[POS_SOURCE];
                               }
 }
 if(i==MIDDLE_POS){
 if   (L[MIDDLE_POS]==1)            derivates_out[POS_TARGET] += -1.*LOSS_MIDDLE*derivates[POS_TARGET];
 if   (R[MIDDLE_POS]==1)            derivates_out[POS_TARGET] += -1.*LOSS_MIDDLE*derivates[POS_TARGET];
 if ( (L[MIDDLE_POS]==0) && (R[MIDDLE_POS]==0)) 
                               {
                                POS_SOURCE = POS_TARGET + POW_2[MIDDLE_POS] + POW_2[N_spins+MIDDLE_POS]; 
                                derivates_out[POS_TARGET] += +2.*LOSS_MIDDLE*derivates[POS_SOURCE];
                               }
 }
// --------------------------------------------------------------------------------------------------
// ------------------------------------- EOF LINDBLAD MIDDLE_POS---------------------------------------------------
// --------------------------------------------------------------------------------------------------

if (TRACE_SWITCH==0)
{
if ( (L[i]==1) && (R[i]==1) ) OUTPUT[97-i] += creal(derivates[POS_TARGET]);

if ( (L[i]==0) && (R[i]==0) ) s_z_temp +=  1.*creal(derivates[POS_TARGET]);
if ( (L[i]==1) && (R[i]==1) ) s_z_temp += -1.*creal(derivates[POS_TARGET]);
OUTPUT[i] += s_z_temp;
if(i%2==0){
	OUTPUT[98] -= (0.5/N_spins)*s_z_temp;
}else{
	OUTPUT[98] += (0.5/N_spins)*s_z_temp;
}
s_z_temp = 0.;
}
if (CURRENT_SWITCH==2){
if( (L[i]==R[i]+1) && (L[i+1]==R[i+1]-1) )  OUTPUT[99]  += NORM*cimag(derivates[POS_TARGET]);
if( (L[i]==R[i]-1) && (L[i+1]==R[i+1]+1) )  OUTPUT[99]  -= NORM*cimag(derivates[POS_TARGET]);
if( (L[i]==R[i]+1) && (L[i+1]==R[i+1]-1) )  OUTPUT[200+i]  += sqrt((L[i+1]+1)*L[i])*cimag(derivates[POS_TARGET]);
if( (L[i]==R[i]-1) && (L[i+1]==R[i+1]+1) )  OUTPUT[200+i]  -= sqrt((L[i]+1)*L[i+1])*cimag(derivates[POS_TARGET]);
//if( (L[1]==1) && (R[1]==0) && (L[2]==0) && (R[2]==1) )  OUTPUT[201]  += 0.5*cimag(derivates[POS_TARGET]);
//if( (L[1]==0) && (R[1]==1) && (L[2]==1) && (R[2]==0) )  OUTPUT[201]  -= 0.5*cimag(derivates[POS_TARGET]);

}
} // end loop Atoms
//printf("end loop atoms\n");
//if (CURRENT_SWITCH==2){
//if( (L[1]==1) && (R[1]==0) && (L[2]==0) && (R[2]==1) )  OUTPUT[201]  += 0.5*cimag(derivates[POS_TARGET]);
//if( (L[1]==0) && (R[1]==1) && (L[2]==1) && (R[2]==0) )  OUTPUT[201]  -= 0.5*cimag(derivates[POS_TARGET]);
//}
  
// ------------------------------------------------------------------------------------------------------
// ------------------------------------- TRACE CHECK ----------------------------------------------------
// ------------------------------------------------------------------------------------------------------
//if (POS_TARGET==0) printf("%.5f\n !!!!!",creal( derivates[POS_TARGET]));
if ( TRACE_SWITCH==0)  {/*printf("k=%i",k);*/ OUTPUT[100]+=creal( derivates[POS_TARGET]);}
} // end loop states

return ;
} // end function derivates

// ###############################################################################
// ###################### Ende Praeambel #########################################  
// ###############################################################################

// ###############################################################################
// #########################  MAIN    ############################################ 
// ###############################################################################

int main(int argc, char* argv[]){  
	
//Zeitmessung
    double time1, timedif;
    time1 = (double) clock();
    time1 = time1 / CLOCKS_PER_SEC;
//Zeitmessung

int i; // site index
int k = 0; // time integration 
double t=0.0;

// ###############################################################################
// ############ CALCULATE THE NUMBER OF DGL - CREATE 2 POWERS #################### 
// ###############################################################################

// here enters the number of states per site!

int *POW_2 = calloc((N_spins+N_spins+2),sizeof(int));
POW_2[0] = 1; 
for (i=1;i<(N_spins+N_spins+2);i++) POW_2[i] = 2*POW_2[i-1];
int N_DGL=POW_2[2*N_spins];
printf("ZDGL=%i -- CALC=%i \n",N_DGL,POW_2[2*N_spins]);

// #################################################################
// ######### START OF INDEX CREATING INDEX ARRAYS ##################
// #################################################################

// use that the density matrix is self adjoint to reduce the numerical effort
// density matrix is written as a vector, so the vector element 
// k =   L[2N_spins-1] * 2^(2N_spins-1) + L[N_spins-2] * 2^(2N_spins-2)  + ... + L[1] * 2^(N_spins+1) + L[0] * 2^N_spins 
//     + R[ N_spins-1] * 2^( N_spins-1) + R[N_spins-2] * 2^( N_spins-2)  + ... + R[1] * 2^1        + R[0] * 2^0  
// of course, the complex conjugate is the vector element k', where L and R Arrays are exchanged.

// #################################################################
// ########### CREATE PARAMETER AND DATA ARRAYS ####################
// #################################################################

complex double *derivates         = calloc( 2*N_DGL,sizeof(double));
complex double *derivates_stst         = calloc( 2*N_DGL,sizeof(double));
complex double *temp              = calloc( 2*N_DGL,sizeof(double));
complex double *k1                = calloc( 2*N_DGL,sizeof(double));
complex double *k2                = calloc( 2*N_DGL,sizeof(double));
// runge kutta vector
// parameter storages
        double *delta_array       = calloc((N_spins+5),sizeof(double));
	    double *coupling_array    = calloc((N_spins+5),sizeof(double));
        double *magnetic_field    = calloc((N_spins+5),sizeof(double));    
	    double *OUTPUT            = calloc( N_OUT,sizeof(double));

#ifdef SHIFT_IMPORT                      	    
	// read in H(i) from MPS-jobfile
	FILE *f1;
	// oeffnen im Lesemodus
	f1 = fopen("shifts.dat", "r");
	if(f1 == NULL) {
		printf("Datei konnte nicht geoeffnet werden.\n");
	}else {
		char weg[200];
		double temp_h;
		fgets(weg,200,f1);
#elif defined(OWN_DISORDER)
	double h_temp = 0.;
	double check_average[N_spins];
	for(int i=0;i<N_spins;i++) check_average[i]=0.;
	srand(time(NULL));
	double current_container=0.0;
	double no_average = 1.*NO_AVERAGE;
	#ifdef H_MAX
	double disorder = H_MAX;
	#else
	double disorder = strtod(argv[1], NULL);
	#endif
	int runno = 0;
for(runno = 1;runno < no_average;runno++){
	printf("start run no %i\n",runno);
#endif
							

for (i=0;i<N_spins;i++)		{
							V(i)     =    V_CONST;
							J(i)     = 1.*J_CONST;
							#ifdef SHIFT_IMPORT
								fscanf(f1,"%lf",&temp_h);
								//if((i==0) || (i==(N_spins-1))){
								H(i) = DISFACTOR*temp_h;
								//}else{
								//H(i) = 2.*DISFACTOR*temp_h;
								//}
								printf("H(%i) = %f\n",i,H(i));
							}
							#elif defined(OWN_DISORDER)
								h_temp = 1.0*rand()/(RAND_MAX); // creates random number in [0.0,1.0) with srand
									//sets h[i] in [-disorder, disorder], as the rng does not cover the negetive part of the interval [-1,1]
									// we need to map [0,1) to the interval [-1,1).
								H(i)  =   ((2.0 * disorder * h_temp) - disorder); 
								printf("disorder at site %i: %.2f\t",i, H(i));
								check_average[i] += H(i); // stores all random values for each site for checking their average values
								printf("check_average at site %i: %.2e\n",i, check_average[i]);
							
							#else
							H(i)      = -1.* H_CONST;
							#endif 
						}
#ifdef SHIFT_IMPORT
fclose(f1);
#endif

for (i=0;i<N_spins;i++) printf("%i: V=%.10f -- J=%.10f -- h=%.10f \n",i,V(i),J(i),H(i));

// #################################################################
// ########### CREATE FILE_NAMES ###################################
// #################################################################
#ifndef OWN_DISORDER
char date[256];time_t curtime;struct tm *loctime;curtime = time (NULL);loctime = localtime (&curtime);
strftime (date, 256, "%m_%d", loctime);
time_t current_time;time(&current_time);struct tm *now;char timestr[1024];now=localtime(&current_time);
strftime(timestr,1024,"%H:%M:%S",now);
char FILE_NAME[2048+1024];
#ifndef PLOT  
snprintf(FILE_NAME,2048+1024,"MEQ_0_1_BASIS_N%i_J%.5f_V%.3f_H%.3f_G_IN%.3f_G_OUT_%.3f.dat",N_spins,J_CONST,V_CONST,H_CONST,PUMP,LOSS); 
#else
snprintf(FILE_NAME,2048+1024,"Dynamikwerte.dat");
#endif
FILE *f_dat;
f_dat=fopen(FILE_NAME,"w");
	fprintf(f_dat,"#t \t ",t); 
	fprintf(f_dat,"N<j> \t <j> \t |<j>|");  
	for(i=0;i<N_spins;i++) { 
	fprintf(f_dat," \t <a+a>%i",i);
	} //falls man an den einzelnen besetzungen interessiert ist
	fprintf(f_dat,"\n"); 
#ifdef CURRENTPLOT
FILE *f_current;
f_current=fopen("current.dat","w");
#endif
#endif

// #################################################################
// ########### END OF CREATE FILE_NAMES ############################
// #################################################################

#ifdef LONGRANGE //setzt Normierung vor longrangeterm
A = Set_A();
#endif

// ###############################################################################
// ############################## INITIAL CONDITIONS #############################
// ###############################################################################
for(int i=0;i<N_DGL;i++) derivates[i]=0.+I*0.;
#if defined(INITIAL_STEADY)
int POS_INIT = 0;
derivates[POS_INIT] = 0.5;
for(int i=0;i<N_spins;i++) POS_INIT += POW_2[i] +POW_2[i+N_spins];
derivates[POS_INIT] = 0.5;
calculate_next_time_step(derivates,k1,t,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
#elif defined(INITIAL_NEEL_0_DOWN)
int POS_INIT = 0;
for(int i=0;i<(int) round(N_spins/2);i++) POS_INIT += POW_2[2*i+1] +POW_2[2*i+1+N_spins];
derivates[POS_INIT] = 1.;
//derivates[0*POW_2[0] + 1*POW_2[1] + 0*POW_2[2] + 1*POW_2[3]+ 0*POW_2[4]+1*POW_2[5]+0*POW_2[6]+1*POW_2[7]] = 1.; // all other states zero 
calculate_next_time_step(derivates,k1,t,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
#elif defined(INITIAL_NEEL_0_UP)
int POS_INIT = 0;
for(int i=0;i<(int) round((N_spins+1)/2);i++) POS_INIT += POW_2[2*i] +POW_2[2*i+N_spins];
derivates[POS_INIT] = 1.;
//derivates[0*POW_2[0] + 1*POW_2[1] + 0*POW_2[2] + 1*POW_2[3]+ 0*POW_2[4]+1*POW_2[5]+0*POW_2[6]+1*POW_2[7]] = 1.; // all other states zero 
calculate_next_time_step(derivates,k1,t,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
#elif defined(INITIAL_ALLUP)
int POS_INIT = 0;
derivates[POS_INIT] = 0.0;
for(int i=0;i<N_spins;i++) POS_INIT += POW_2[i] +POW_2[i+N_spins];
derivates[POS_INIT] = 1.0;
calculate_next_time_step(derivates,k1,t,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
#elif defined(ALL_DOWN_LAST_UP)
int POS_INIT = POW_2[N_spins-1]+POW_2[N_spins+N_spins-1];
derivates[POS_INIT] = 1.0;
#else
derivates[0]=1.;
#endif
for(int i=0;i<1000;i++) OUTPUT[i]=0.;
OUTPUT[100] =1.;
OUTPUT[999] =10.;
// ###############################################################################
// ####################### END OF INITIAL CONDITIONS #############################
// ###############################################################################

// ###########################################################################
// ########################## TIME _DOMAIN SOLUTION ##########################
// ###########################################################################
//Beginn der Integration von 0fs bis N*Schrittweite in fs
	k=0;
   while (k<N && OUTPUT[999]>MAX_DEVIA)
     {  
       t=delta_t*k;
       // #######################################
      if  (k % STEP_IS_SAVED ==0) 
         {
		  #ifndef OWN_DISORDER
			  fprintf(f_dat,"%e \t ",t); 
			  fprintf(f_dat,"%e \t %e \t %e",OUTPUT[98],N_spins*OUTPUT[99],cabs(OUTPUT[98]));  
			  for(i=0;i<N_spins;i++) { 
				fprintf(f_dat," \t %e",creal(OUTPUT[97-i]));
			  } //falls man an den einzelnen sz interessiert ist
			  #ifdef CURRENTPLOT
			  fprintf(f_current,"%e",t);
			  for(i=0;i<(N_spins);i++) fprintf(f_current,"\t %e",OUTPUT[200+i]);
			  fprintf(f_current,"\n");
			  #endif
			  fprintf(f_dat,"\n"); 
          #endif
         }
      if (cabs(OUTPUT[100]-1.000)<0.10) 
      {
// ################################################################
// ############### 4 VECTOR RUNGE KUTTA ###########################
// ################################################################

calculate_next_time_step(derivates,k1,t,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
CALCULATE_TEMP(derivates,k1,delta_t*HALFTE,temp,i);

calculate_next_time_step(temp, k2, t+delta_t*0.5,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
CALCULATE_TEMP(derivates,k2,delta_t*HALFTE,temp,i);ADD_VECTOR(k1,k2,2.0,i);

calculate_next_time_step(temp, k2, t+delta_t*0.5,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
CALCULATE_TEMP(derivates,k2,delta_t,temp,i);ADD_VECTOR(k1,k2,2.0,i);

calculate_next_time_step(temp, k2, t+delta_t,delta_array,coupling_array,magnetic_field,POW_2,OUTPUT); 
ADD_VECTOR(k1,k2,1.0,i);ADD_VECTOR(derivates,k1,delta_t*SECHSTEL,i);

// ################################################################
// ############ END OF 4 VECTOR RUNGE KUTTA #######################
// ################################################################
     } else {
		printf("trace zu weit weg\n");break; 
	 }
     STEADY_STATE_TEST(derivates,derivates_stst,OUTPUT,N_DGL); //Berechne Abstand zum StSt
     #ifndef OWN_DISORDER
     if (k % TIME_OUTPUT == 0){ 
		printf("Jt=%.1f - stst_abw %.10f - TR=%.5f",t,OUTPUT[999],OUTPUT[100]);
		for(i=0;i<N_spins;i++) { 
		printf(" -- Site[%i]_occ=%.5f",i,creal(OUTPUT[97-i]));
		}
		printf("-- <M> %.5f -- <j> = %.5f-- N<j> = %.5f-- <j>_0 = %.5f\n",OUTPUT[98],OUTPUT[99],N_spins*OUTPUT[99],N_spins*OUTPUT[201]);
	}
	#endif
	k++;
}// ende while
        
  printf("CHECK FOR 1 SPIN: SteadyState ANA Sz=%.10f -- NUM Sz=%.10f \n",(LOSS-PUMP)/(LOSS+PUMP),OUTPUT[0]);       
  #ifndef OWN_DISORDER
  fclose(f_dat);
	#ifdef PLOT
	system("gnuplot -p plot.gp");	
	#endif
  #ifdef CURRENTPLOT
  fclose(f_current);
  #endif
  #else
	current_container += OUTPUT[99];
  	}//end of loop for runno
	current_container;
	printf("average absolute current after %i runs is %f, with h_max = %.2f and an average disorder of\n",runno,current_container/no_average,disorder);
	for(int i=0;i<N_spins;i++) printf("site %i: %f\n",i,check_average[i]/no_average);
	FILE *f_disorder;
	char FILE_NAME[2048+1024];
	snprintf(FILE_NAME,2048+1024,"RK_%iruns_hmax_%.2f",runno,disorder); 
	f_disorder=fopen(FILE_NAME,"w");
	fprintf(f_disorder,"average absolute current after %i runs is %f, with h_max = %.2f and an average disorder of\n",runno,current_container/no_average,disorder);
	for(int i=0;i<N_spins;i++) fprintf(f_disorder,"site %i: %f\n",i,check_average[i]/no_average);
	fclose(f_disorder);
  #endif



  		  
  free(coupling_array);free(delta_array);free(magnetic_field);
  free(temp);free(derivates);free(k1);free(k2);
  free(POW_2);free(OUTPUT);

  //Zeitmessung
    timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - time1;
    printf("Runtime: %.2f minutes\n", timedif/60.0);
	//Zeitmessung


  return 0;
}

