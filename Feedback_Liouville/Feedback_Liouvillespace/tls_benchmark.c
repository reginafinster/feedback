 #include <math.h>
 #include <stdlib.h>
 #include <complex.h>
 #include <stdio.h>
 #include <time.h>
 
//PARAMETERBEREICH

#define N_steps 1000  	 	//Anzahl Zeitschritte
#define t0 0.0				//Zeit - Anfang
#define t1 100.0			//Ende


#define INIT_EXCITATION 1.0 // Für künstl. Initial-state-Anregung. für 0.0 alles im GZS, 1.0 alles im ex.st.


// WARNING ich muss die Gammas Wurzelziehen im vgl zu MPS Programm 
#define gamma_in 0.0 	//Lindblad-Pump
#define gamma_out 0.5	//Lindblad-Loss
#define kappa 1.0		//Dephasing
#define omega 0.0		//Coherent Driving

//PARAMETERBEREICH ENDE

#define Anzahl_Dgl (2*2) //je 2 pro ZNS und h.c.

//Überschreibt in derivates den Eintrag, der für diese Kombination der Zustände im Speicher liegt
//(a,b = bra,ket von TLS)
#define rho(a,b) ( *(derivates  + (a)*2 + (b) ))
#define rho_out(a,b) ( *(derivates_out  + (a)*2 + (b) ))


// RK-Makros
#define calc_temp(tmp, a, b, fak, l) ({for(l=0;l<N_dgl;l++) (tmp)[l]=(a)[l]+(fak)*(b)[l];}) 	// xi + fak*ki
#define add_vect(tmp, neu, fak, l) ({for(l=0;l<N_dgl;l++) (tmp)[l]+=(fak)*(neu)[l];})			// vektoradd (für sechstel*(k0+2k1...
#define sechstel 0.1666666666666666666666667
#define haelfte	0.5

void calc_step (
complex double *derivates, //array aus main
complex double *derivates_out,  //array an main zurück
        double t,
        double *OUTPUT,
        int k
)
{
for (int x=0;x<50;x++) OUTPUT[x]=0.; //OUTPUT - Reset

for (int a=0;a<2;a++)
{
for (int b=0;b<2;b++)
{ 

rho_out(a,b) = 0. + I * 0.;  //rho_out - Reset

//MASTERGLEICHUNG//

//Lindblad Loss
							rho_out(a,b) 	+= -0.5*gamma_out*(a+b)*rho(a,b);
if ( (a==0) && (b==0)  )	rho_out(a,b) 	+= gamma_out*rho(a+1,b+1);

//Lindblad Pump
if ( (a==0) )				rho_out(a,b) 	+= -0.5*gamma_in*rho(a,b);
if ( (b==0) )				rho_out(a,b) 	+= -0.5*gamma_in*rho(a,b);
if ( (a==1) && (b==1)  )	rho_out(a,b) 	+= gamma_in*rho(a-1,b-1);

//Dephasing
if ( (a==1) && (b==0)  )	rho_out(a,b) 	+= -0.5*kappa*rho(a,b); //man könnte auch einfach fordern "a!=b"...
if ( (a==0) && (b==1)  )	rho_out(a,b) 	+= -0.5*kappa*rho(a,b);

//Coherent Driving
if ( (a==1) )				rho_out(a,b)	+= -I*omega*rho(a-1,b);
if ( (a==0) )				rho_out(a,b)	+= -I*omega*rho(a+1,b);
if ( (b==1) )				rho_out(a,b)	+= I*omega*rho(a,b-1);
if ( (b==0) )				rho_out(a,b)	+= I*omega*rho(a,b+1);

//MASTERGLEICHUNG//

if ((a==b)){
	OUTPUT[0] += creal(rho(a,b)) ; 			// Spur - sollte 1 ergeben
	if (a==0) OUTPUT[1] += creal(rho(a,b)); 	// Besetzung Grundzustand
	if (a==1) OUTPUT[2] += creal(rho(a,b)); 	// Besetzung angeregter Zustand
}

}//ende for loop b
}//ende for loop a
return;
}//ende calc-step function

int main()
{
	
//Zeitmessung
    double time1, timedif;
    time1 = (double) clock();
    time1 = time1 / CLOCKS_PER_SEC;
//Zeitmessung
	
double delta_t = (t1-t0)/N_steps;			//Zeitschritt - "h"

int N_dgl = Anzahl_Dgl; 					//Anzahl DGL 
	
printf("delta_t=%f\n",delta_t);
printf("N_dgl=%i\n",N_dgl);

	
//Speicher reservieren	
	
complex double *derivates = calloc(2*N_dgl,sizeof(double));
// Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
complex double *temp = calloc(2*N_dgl,sizeof(double));
// Die errechneten Zwischenfunktionswerte - die Steigungen
complex double *k1     = calloc(2*N_dgl,sizeof(double));
complex double *k2     = calloc(2*N_dgl,sizeof(double));
        double *OUTPUT = calloc(100,sizeof(double)); 
     
FILE *f_dat;
f_dat=fopen("Funktionswerte.dat","w");

for (int x=0;x<N_dgl;x++) derivates[x] = 0. + I* 0.;    //Setze derivates auf 0

OUTPUT[0] = 1.;										//Setze Spur auf 1

//Anfangswahrscheinlichkeiten: Summe = 1, rho(a,b)
rho(1,1)= INIT_EXCITATION;
rho(0,0)= 1.-INIT_EXCITATION;


int l; // für RK-Makros
double t; // Zeit

//RK 4 - LOOP//
for (int k=0; k<(N_steps+1); k++) {
t=delta_t*k;	//Zeit -> 'Abszisse'
	
calc_step(derivates, k1, t, OUTPUT, k); 		//1. Schritt
calc_temp(temp, derivates, k1, delta_t*haelfte,l); 		//1. Steigung berechnen + h/2

calc_step(temp, k2, t+delta_t*0.5,OUTPUT,k);	//2. Schritt
calc_temp(temp, derivates, k2, delta_t*haelfte,l);		//2. Steigung berechnen + h/2
add_vect(k1,k2,2.0,l);									//k2 zweimal auf k1

calc_step(temp, k2, t+delta_t*0.5,OUTPUT,k);	//3. Schritt
calc_temp(temp, derivates, k2, delta_t,l);		    	//3. Steigung berechnen + h
add_vect(k1,k2,2.0,l);									//k2 zweimal auf k1

calc_step(temp, k2, t+delta_t,OUTPUT,k);		//4. Schritt
add_vect(k1,k2,1.0,l);									//k2 einmal auf k1
add_vect(derivates,k1,delta_t*sechstel,l);				//1/6 k1 auf derivates
if(k % (N_steps/20) == 0){
printf("%0.0f Prozent fertig. -- TRACE = %.10f -- rho11 = %.5f -- rho00 = %.5f \n",k*100./N_steps,OUTPUT[0],OUTPUT[2],OUTPUT[1]);
}
fprintf(f_dat,"%.10f \t %.10f \t %.10f \n",t,OUTPUT[1],OUTPUT[2]);
}
//RK 4 - LOOP//


// Speicher freigeben

  fclose(f_dat);
  free(temp);
  free(derivates);
  free(k1);
  free(k2);
  free(OUTPUT);
  
//Zeitmessung
    timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - time1;
    printf("Finished. Runtime: %.2f minutes\n", timedif/60.0);
//Zeitmessung

  return 0;
}
