

#include "QST_parameter.h"


// Header for Calculation of Eigenvalues
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


// Zur Definition der Zustaende
#define N_l (6) // Zahl der Level
#define rho(a,b)      ( *(derivates     + (a)*N_l + (b)  ))  // Werte der Einträge der Dichtematrix
#define rho_out(a,b)  ( *(derivates_out + (a)*N_l + (b)  ))  // Ableitung der Dichtematrix
//#define rho_stst(a,b) ( *(derivates_stst+ (a)*N_l + (b)  ))


#define Zahl_der_DGL      ( N_l*N_l + 10  )
int N_DGL=Zahl_der_DGL;


#define PURE_DEPHASING
#define CROSS_DEPHASING
#define KAPPA ( 0.050 )

//#define CSPINFLIP
//#define ISPINFLIP



// -----------------------------------------------------------------------
// -----  Analytical joint photodetection probability (no spinflip)  -----
// -----------------------------------------------------------------------
double P_ana(double theta1, double phi1, double theta2, double phi2, double FSS, double t)
{

    double p;
    // t*=1000;

    p = 0.5 * exp(-GAM_R * t) * cabs( cos(0.5 * (theta1-theta2)) * cos( 0.5 *(phi1+phi2 +1.0*FSS*t)) + I * cos( 0.5 * (theta1+theta2)) * sin( 0.5 * (phi1+phi2 +1.0*FSS*t)) ) * cabs( cos(0.5 * (theta1-theta2)) * cos( 0.5 *(phi1+phi2 +1.0*FSS*t)) + I * cos( 0.5 * (theta1+theta2)) * sin( 0.5 * (phi1+phi2 +1.0*FSS*t)) );

    return p;
}


// #################################################################################
// ############## DGL system (Mastergleichung)  ######################################################
// #################################################################################

void calculate_next_time_step
(
 complex double *derivates,
 complex double *derivates_out,
 double t,
 double *OUTPUT,
 double FSS,
 double FLIP
 )
{


    int a,b;

    for (a=0;a<100;a++) OUTPUT[a]=0.;


    for (a=0;a<N_l;a++)
    {
        for (b=0;b<N_l;b++)
        {
            rho_out(a,b) = 0. + I*0.;  // Set derivatives back to zero

            // ----------------------------------------------------------------------
            // ------------- Left part - Density Matrix -----------------------------
            // ----------------------------------------------------------------------
            // ##### Ground state
            if  (a==0)            {rho_out(a,b) += 0. * rho(a,b) ;}

            // ##### V-biexciton photon x V-exciton state
            if  (a==1)            {rho_out(a,b) += -0.5 * GAM_R_XV                     * rho(a,b)
                                                   -I * 0.5 * FSS                      * rho(a,b)  ;}

            // #####  V-biexciton photon x H-exciton state
            if  (a==2)            {rho_out(a,b) += -0.5 * GAM_R_XH                     * rho(a,b)
                                                   +I * 0.5 * FSS                      * rho(a,b)  ;}

            // ##### H-biexciton photon x V-exciton state
            if  (a==3)            {rho_out(a,b) += -0.5 * GAM_R_XV                     * rho(a,b)
                                                   -I * 0.5 * FSS                      * rho(a,b)  ;}

            // ##### H-biexciton photon x H-exciton state
            if  (a==4)            {rho_out(a,b) += - 0.5 * GAM_R_XH                    * rho(a,b)
                                                   +I * 0.5 * FSS                      * rho(a,b)  ;}

            // ##### Biexciton state
            if  (a==5)            {rho_out(a,b) += - 0.5 * (GAM_R_BH + GAM_R_BV)       * rho(a,b)
                                                   +I * 0.5 * DELTA_XX                 * rho(a,b)  ;}

            // ----------------------------------------------------------------------
            // ------------- Right part - Density Matrix -----------------------------
            // ----------------------------------------------------------------------
            // ##### ground state
            if  (b==0)            {rho_out(a,b) += 0. * rho(a,b) ;}

            // ##### V-biexciton photon x V-exciton state
            if  (b==1)            {rho_out(a,b) += -0.5 * GAM_R_XV                      * rho(a,b)
                                                   +I * 0.5 * FSS                       * rho(a,b)  ;}

            // ##### V-biexciton photon x H-exciton state
            if  (b==2)            {rho_out(a,b) += -0.5 * GAM_R_XH                      * rho(a,b)
                                                   -I * 0.5 * FSS                       * rho(a,b)  ;}

            // ##### H-biexciton photon x V-exciton state
            if  (b==3)            {rho_out(a,b) += -0.5 * GAM_R_XV                      * rho(a,b)
                                                   +I * 0.5 * FSS                       * rho(a,b)  ;}

            // ##### H-biexciton photon x H-exciton state
            if  (b==4)            {rho_out(a,b) += -0.5 * GAM_R_XH                      * rho(a,b)
                                                   -I * 0.5 * FSS                       * rho(a,b)  ;}

            // ##### Biexciton state
            if  (b==5)            {rho_out(a,b) += -0.5 * (GAM_R_BH + GAM_R_BV)         * rho(a,b)
                                                   -I * 0.5 * DELTA_XX                  * rho(a,b)  ;}
            // ----------------------------------------------------------------------
            // -----------------  Lindblad Radiative Recombination  -----------------
            // ----------------------------------------------------------------------
            if ((a==0)&&(b==0)) rho_out(a,b) += +1.* GAM_R_XV * rho(1,1)
                                                +1.* GAM_R_XH * rho(2,2)
                                                +1.* GAM_R_XV * rho(3,3)
                                                +1.* GAM_R_XH * rho(4,4);

            /*
            if ((a==1)&&(b==1)) rho_out(a,b) += +0.5* GAM_R_BV * rho(5,5);
            if ((a==2)&&(b==2)) rho_out(a,b) += +0.5* GAM_R_BH * rho(5,5);
            if ((a==3)&&(b==3)) rho_out(a,b) += +0.5* GAM_R_BV * rho(5,5);
            if ((a==4)&&(b==4)) rho_out(a,b) += +0.5* GAM_R_BH * rho(5,5);
           */


            // ----------------------------------------------------------------------
            // ---------------  Lindblad Pure dephasing  ----------------------------
            // ----------------------------------------------------------------------
#ifdef PURE_DEPHASING
            if ((a==1)&&(b==0)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==0)&&(b==1)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==2)&&(b==0)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);
            if ((a==0)&&(b==2)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);
            if ((a==3)&&(b==0)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==0)&&(b==3)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==4)&&(b==0)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);
            if ((a==0)&&(b==4)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);

            if ((a==1)&&(b==5)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==5)&&(b==1)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==2)&&(b==5)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);
            if ((a==5)&&(b==2)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);
            if ((a==3)&&(b==5)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==5)&&(b==3)) rho_out(a,b) += -1.* GAM_PV * rho(a,b);
            if ((a==4)&&(b==5)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);
            if ((a==5)&&(b==4)) rho_out(a,b) += -1.* GAM_PH * rho(a,b);

#ifdef CROSS_DEPHASING
            if ((a==3)&&(b==4)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
            if ((a==4)&&(b==3)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
            if ((a==3)&&(b==2)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
            if ((a==2)&&(b==3)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
            if ((a==1)&&(b==4)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
            if ((a==4)&&(b==1)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
            if ((a==1)&&(b==2)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
            if ((a==2)&&(b==1)) rho_out(a,b) += -0.5 * KAPPA * (GAM_PH + GAM_PV) * rho(a,b);
#endif

#endif


            // ----------------------------------------------------------------------
            // ----------------  Coherent Polarization flip  ------------------------
            // ----------------------------------------------------------------------
#ifdef CSPINFLIP
            if (a==1)            {rho_out(a,b) += -I * FLIP * rho(2,b);}
            if (a==2)            {rho_out(a,b) += -I * FLIP * rho(1,b);}
            if (a==3)            {rho_out(a,b) += -I * FLIP * rho(4,b);}
            if (a==4)            {rho_out(a,b) += -I * FLIP * rho(3,b);}

            if (b==1)            {rho_out(a,b) += +I * FLIP * rho(a,2);}
            if (b==2)            {rho_out(a,b) += +I * FLIP * rho(a,1);}
            if (b==3)            {rho_out(a,b) += +I * FLIP * rho(a,4);}
            if (b==4)            {rho_out(a,b) += +I * FLIP * rho(a,3);}
#endif

            // ----------------------------------------------------------------------
            // ---------------  Incoherent Polarization flip  -----------------------
            // ----------------------------------------------------------------------
#ifdef ISPINFLIP
            if (a==1)            {rho_out(a,b) += -0.5 * FLIP * rho(1,b);}
            if (a==2)            {rho_out(a,b) += -0.5 * FLIP * rho(2,b);}
            if (a==3)            {rho_out(a,b) += -0.5 * FLIP * rho(3,b);}
            if (a==4)            {rho_out(a,b) += -0.5 * FLIP * rho(4,b);}

            if (b==1)            {rho_out(a,b) += -0.5 * FLIP * rho(a,1);}
            if (b==2)            {rho_out(a,b) += -0.5 * FLIP * rho(a,2);}
            if (b==3)            {rho_out(a,b) += -0.5 * FLIP * rho(a,3);}
            if (b==4)            {rho_out(a,b) += -0.5 * FLIP * rho(a,4);}

            if( (a==1) && (b==1) )  {rho_out(a,b) += +1.0 * FLIP * rho(2,2);}
            if( (a==2) && (b==2) )  {rho_out(a,b) += +1.0 * FLIP * rho(1,1);}
            if( (a==3) && (b==3) )  {rho_out(a,b) += +1.0 * FLIP * rho(4,4);}
            if( (a==4) && (b==4) )  {rho_out(a,b) += +1.0 * FLIP * rho(3,3);}
#endif


            // ----------------------------------------------------------------------
            // -------------  Expectation Values  -----------------------------------
            // ----------------------------------------------------------------------
            if ( (a==b)              ) OUTPUT[0] += creal(rho(a,b));


            if ( (a==b)  && (a==0)   ) OUTPUT[6] += creal(rho(a,b));
            if ( (a==b)  && (a==1)   ) OUTPUT[1] += creal(rho(a,b));
            if ( (a==b)  && (a==2)   ) OUTPUT[2] += creal(rho(a,b));
            if ( (a==b)  && (a==3)   ) OUTPUT[3] += creal(rho(a,b));
            if ( (a==b)  && (a==4)   ) OUTPUT[4] += creal(rho(a,b));
            if ( (a==b)  && (a==5)   ) OUTPUT[5] += creal(rho(a,b));

        } //b
    } //a

    return ;
}

// ###############################################################################
// ###############################################################################
// #########################  MAIN    ############################################
// ###############################################################################
// ###############################################################################
int main()
{
    int a,b,k,n,i,z;
    double t;

/*
#ifdef ENTANGLEMENT
    char FILE_NAME[2048+1024];
    snprintf(FILE_NAME,2048+1024,"entanglement_PUMP%.5f.txt",1000.*PUMP);
    FILE *f_ent;
    f_ent=fopen(FILE_NAME,"w");
#endif
*/

    // #############################################################################
    // #############################  PREPARATION RK4  #############################
    // #############################################################################
    complex double *derivates    = calloc(2*N_DGL,sizeof(double));
    // Um den Anfangswert zu speichern, mit dem die Zwischen-Anfangswerte berechnet werden
    complex double *temp = calloc(2*N_DGL,sizeof(double));
    // Die errechneten Zwischenfunktionswerte - die Steigungen
    complex double *k1 = calloc(2*N_DGL,sizeof(double));
    complex double *k2 = calloc(2*N_DGL,sizeof(double));
    double *OUTPUT = calloc(100,sizeof(double));




        // Create Files for time traces of the 36 joint photodetection probabilities (fold = convolved with detector function, ana = analytical result without spinflip)

    FILE *f_pro_HX;
    f_pro_HX = fopen("pro_HX.txt","w");
    FILE *f_pro_VX;
    f_pro_VX = fopen("pro_VX.txt","w");
    FILE *f_pro_DX;
    f_pro_DX = fopen("pro_DX.txt","w");
    FILE *f_pro_AX;
    f_pro_AX = fopen("pro_AX.txt","w");
    FILE *f_pro_RX;
    f_pro_RX = fopen("pro_RX.txt","w");
    FILE *f_pro_LX;
    f_pro_LX = fopen("pro_LX.txt","w");

    FILE *f_pro_fold_HX;
    f_pro_fold_HX = fopen("pro_fold_HX.txt","w");
    FILE *f_pro_fold_VX;
    f_pro_fold_VX = fopen("pro_fold_VX.txt","w");
    FILE *f_pro_fold_DX;
    f_pro_fold_DX = fopen("pro_fold_DX.txt","w");
    FILE *f_pro_fold_AX;
    f_pro_fold_AX = fopen("pro_fold_AX.txt","w");
    FILE *f_pro_fold_RX;
    f_pro_fold_RX = fopen("pro_fold_RX.txt","w");
    FILE *f_pro_fold_LX;
    f_pro_fold_LX = fopen("pro_fold_LX.txt","w");
    
    FILE *fit_HH;
    fit_HH = fopen("fit_HH.txt", "w");
    FILE *fit_HV;
    fit_HV = fopen("fit_HV.txt", "w");
    FILE *fit_HD;
    fit_HD = fopen("fit_HD.txt", "w");
    FILE *fit_HR;
    fit_HR = fopen("fit_HR.txt", "w");
    FILE *fit_VH;
    fit_VH = fopen("fit_VH.txt", "w");
    FILE *fit_VV;
    fit_VV = fopen("fit_VV.txt", "w");
    FILE *fit_VD;
    fit_VD = fopen("fit_VD.txt", "w");
    FILE *fit_VR;
    fit_VR = fopen("fit_VR.txt", "w");
    FILE *fit_DH;
    fit_DH = fopen("fit_DH.txt", "w");
    FILE *fit_DV;
    fit_DV = fopen("fit_DV.txt", "w");
    FILE *fit_DD;
    fit_DD = fopen("fit_DD.txt", "w");
    FILE *fit_DR;
    fit_DR = fopen("fit_DR.txt", "w");
    FILE *fit_RH;
    fit_RH = fopen("fit_RH.txt", "w");
    FILE *fit_RV;
    fit_RV = fopen("fit_RV.txt", "w");
    FILE *fit_RD;
    fit_RD = fopen("fit_RD.txt", "w");
    FILE *fit_RR;
    fit_RR = fopen("fit_RR.txt", "w");
    
    
    FILE *f_pro_ana_HX;
    f_pro_ana_HX = fopen("pro_ana_HX.txt","w");
    FILE *f_pro_ana_VX;
    f_pro_ana_VX = fopen("pro_ana_VX.txt","w");
    FILE *f_pro_ana_DX;
    f_pro_ana_DX = fopen("pro_ana_DX.txt","w");
    FILE *f_pro_ana_AX;
    f_pro_ana_AX = fopen("pro_ana_AX.txt","w");
    FILE *f_pro_ana_RX;
    f_pro_ana_RX = fopen("pro_ana_RX.txt","w");
    FILE *f_pro_ana_LX;
    f_pro_ana_LX = fopen("pro_ana_LX.txt","w");



    FILE *f_neg;
    FILE *f_av_neg;
    f_neg = fopen("negativity.txt","w");
    f_av_neg = fopen("av_neg.txt", "w");

    // Definition of System Parameters
    // Here, you can add for/while loops in case you want to study a certain parameter range

    double FLIP;
    FLIP = 0.011;

    double FSS;
    FSS = + (0.0128 +0.0000) / hbar ;


    // time-resolved probabilities

    double hh[MAX_STEPS],hv[MAX_STEPS],hd[MAX_STEPS],ha[MAX_STEPS],hr[MAX_STEPS],hl[MAX_STEPS];
    double vh[MAX_STEPS],vv[MAX_STEPS],vd[MAX_STEPS],va[MAX_STEPS],vr[MAX_STEPS],vl[MAX_STEPS];
    double dh[MAX_STEPS],dv[MAX_STEPS],dd[MAX_STEPS],da[MAX_STEPS],dr[MAX_STEPS],dl[MAX_STEPS];
    double ah[MAX_STEPS],av[MAX_STEPS],ad[MAX_STEPS],aa[MAX_STEPS],ar[MAX_STEPS],al[MAX_STEPS];
    double rh[MAX_STEPS],rv[MAX_STEPS],rd[MAX_STEPS],ra[MAX_STEPS],rr[MAX_STEPS],rl[MAX_STEPS];
    double lh[MAX_STEPS],lv[MAX_STEPS],ld[MAX_STEPS],la[MAX_STEPS],lr[MAX_STEPS],ll[MAX_STEPS];




    // helpers to calculate negativity

    double PTrho[32];  //partially transposed matrix
    double value[4] = { 0. };
    double sum;

    double NEG[MAX_STEPS];

    gsl_matrix_complex_view m;
    gsl_vector *eval = gsl_vector_calloc(4);
    gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(4);

    // negativity averaged over time resolution of experimental setup

    double avePTrho[32][MAX_STEPS] = { 0. }; // 100/delta_t = 200
    double ave_value[4][MAX_STEPS] = { 0. };
    double ave_NEG[MAX_STEPS]={0.};

    // ###############################################################################
    // #############################  INITIAL CONDITION  #############################
    // ###############################################################################

    for (n=0;n<N_DGL;n++) {
        derivates[n]=0.;
    }

    for(a=0;a<N_l;a++){
        for(b=0;b<N_l;b++){
        // |Psi(t=t_XX=0)> = (|HH> +|VV>)/Sqrt(2)
        // -----------
        if ((a==4)&&(b==4)) rho(a,b) = 0.5;
        if ((a==1)&&(b==1)) rho(a,b) = 0.5;
        // ---------
        if ((a==4)&&(b==1)) rho(a,b) = 0.0;
        if ((a==1)&&(b==4)) rho(a,b) = 0.0;
        // -----------
    }}







    k=0; // initialize k
    // ########################## TAU_DOMAIN SOLUTION ############################

    //printf("t \t HB \t HV \t HH \t HG \t g2HH \t g2HV \t VB \t VV \t VH \t VG \t g2VH \t g2VV \n");
    //printf("t \t g2HH \t g2HV \t g2VH \t g2VV \n");

    // integrate in tau until steady state
    while(k<MAX_STEPS-1)
    {
        k++;
        t=delta_t*k;


        // ############### 4 VECTOR RUNGE KUTTA ###########################
        calculate_next_time_step(derivates, k1 , t,OUTPUT,FSS,FLIP);
        CALCULATE_TEMP(derivates,k1,delta_t*HALFTE,temp,i);
        calculate_next_time_step(temp, k2, t+delta_t*0.5,OUTPUT,FSS,FLIP);
        CALCULATE_TEMP(derivates,k2,delta_t*HALFTE,temp,i);
        ADD_VECTOR(k1,k2,2.0,i);
        calculate_next_time_step(temp, k2, t+delta_t*0.5,OUTPUT,FSS,FLIP);
        CALCULATE_TEMP(derivates,k2,delta_t,temp,i);
        ADD_VECTOR(k1,k2,2.0,i);
        calculate_next_time_step(temp, k2, t+delta_t,OUTPUT,FSS,FLIP);
        ADD_VECTOR(k1,k2,1.0,i);
        ADD_VECTOR(derivates,k1,delta_t*SECHSTEL,i);


        hh[k] = 1.00 * creal(rho(4,4));
        hv[k] = 1.00 * creal(rho(3,3));
        hd[k] = 0.50 * creal(rho(4,4) +rho(4,3) +rho(3,4) +rho(3,3));
        ha[k] = 0.50 * creal(rho(4,4) -rho(4,3) -rho(3,4) +rho(3,3));
        hr[k] = 0.50 * creal(rho(4,4) +I*rho(4,3) -I*rho(3,4) +rho(3,3));
        hl[k] = 0.50 * creal(rho(4,4) -I*rho(4,3) +I*rho(3,4) +rho(3,3));

        vh[k] = 1.00 * creal(rho(2,2));
        vv[k] = 1.00 * creal(rho(1,1));
        vd[k] = 0.50 * creal(rho(2,2) +rho(2,1) +rho(1,2) +rho(1,1));
        va[k] = 0.50 * creal(rho(2,2) -rho(2,1) -rho(1,2) +rho(1,1));
        vr[k] = 0.50 * creal(rho(2,2) +I*rho(2,1) -I*rho(1,2) +rho(1,1));
        vl[k] = 0.50 * creal(rho(2,2) -I*rho(2,1) +I*rho(1,2) +rho(1,1));

        dh[k] = 0.50 * creal(rho(4,4) +rho(4,2) +rho(2,4) +rho(2,2));
        dv[k] = 0.50 * creal(rho(3,3) +rho(3,1) +rho(1,3) +rho(1,1));
        dd[k] = 0.25 * creal( rho(4,4) +rho(4,3) +rho(4,2) +rho(4,1)
                             +rho(3,4) +rho(3,3) +rho(3,2) +rho(3,1)
                             +rho(2,4) +rho(2,3) +rho(2,2) +rho(2,1)
                             +rho(1,4) +rho(1,3) +rho(1,2) +rho(1,1));
        da[k] = 0.25 * creal( rho(4,4) -rho(4,3) +rho(4,2) -rho(4,1)
                             -rho(3,4) +rho(3,3) -rho(3,2) +rho(3,1)
                             +rho(2,4) -rho(2,3) +rho(2,2) -rho(2,1)
                             -rho(1,4) +rho(1,3) -rho(1,2) +rho(1,1));
        dr[k] = 0.25 * creal( rho(4,4) +I*rho(4,3) +rho(4,2) +I*rho(4,1)
                             -I*rho(3,4) +rho(3,3) -I*rho(3,2) +rho(3,1)
                             +rho(2,4) +I*rho(2,3) +rho(2,2) +I*rho(2,1)
                             -I*rho(1,4) +rho(1,3) -I*rho(1,2) +rho(1,1));
        dl[k] = 0.25 * creal( rho(4,4) -I*rho(4,3) +rho(4,2) -I*rho(4,1)
                             +I*rho(3,4) +rho(3,3) +I*rho(3,2) +rho(3,1)
                             +rho(2,4) -I*rho(2,3) +rho(2,2) -I*rho(2,1)
                             +I*rho(1,4) +rho(1,3) +I*rho(1,2) +rho(1,1));

        ah[k] = 0.50 * creal(rho(4,4) -rho(4,2) -rho(2,4) +rho(2,2));
        av[k] = 0.50 * creal(rho(3,3) -rho(3,1) -rho(1,3) +rho(1,1));
        ad[k] = 0.25 * creal( rho(4,4) +rho(4,3) -rho(4,2) -rho(4,1)
                             +rho(3,4) +rho(3,3) -rho(3,2) -rho(3,1)
                             -rho(2,4) -rho(2,3) +rho(2,2) +rho(2,1)
                             -rho(1,4) -rho(1,3) +rho(1,2) +rho(1,1));
        aa[k] = 0.25 * creal( rho(4,4) -rho(4,3) -rho(4,2) +rho(4,1)
                             -rho(3,4) +rho(3,3) +rho(3,2) -rho(3,1)
                             -rho(2,4) +rho(2,3) +rho(2,2) -rho(2,1)
                             +rho(1,4) -rho(1,3) -rho(1,2) +rho(1,1));
        ar[k] = 0.25 * creal( rho(4,4) +I*rho(4,3) -rho(4,2) -I*rho(4,1)
                             -I*rho(3,4) +rho(3,3) +I*rho(3,2) -rho(3,1)
                             -rho(2,4) -I*rho(2,3) +rho(2,2) +I*rho(2,1)
                             +I*rho(1,4) -rho(1,3) -I*rho(1,2) +rho(1,1));
        al[k] = 0.25 * creal( rho(4,4) -I*rho(4,3) -rho(4,2) +I*rho(4,1)
                             +I*rho(3,4) +rho(3,3) -I*rho(3,2) -rho(3,1)
                             -rho(2,4) +I*rho(2,3) +rho(2,2) -I*rho(2,1)
                             -I*rho(1,4) -rho(1,3) +I*rho(1,2) +rho(1,1));

        rh[k] = 0.50 * creal(rho(4,4) +I*rho(4,2) -I*rho(2,4) +rho(2,2));
        rv[k] = 0.50 * creal(rho(3,3) +I*rho(3,1) -I*rho(1,3) +rho(1,1));

        rd[k] = 0.25 * creal( rho(4,4) +rho(4,3) +I*rho(4,2) +I*rho(4,1)
                             +rho(3,4) +rho(3,3) +I*rho(3,2) +I*rho(3,1)
                             -I*rho(2,4) -I*rho(2,3) +rho(2,2) +rho(2,1)
                             -I*rho(1,4) -I*rho(1,3) +rho(1,2) +rho(1,1));

        ra[k] = 0.25 * creal( rho(4,4) -rho(4,3) +I*rho(4,2) -I*rho(4,1)
                             -rho(3,4) +rho(3,3) -I*rho(3,2) +I*rho(3,1)
                             -I*rho(2,4) +I*rho(2,3) +rho(2,2) -rho(2,1)
                             +I*rho(1,4) -I*rho(1,3) -rho(1,2) +rho(1,1));

        rr[k] = 0.25 * creal( rho(4,4) +I*rho(4,3) +I*rho(4,2) -rho(4,1)
                             -I*rho(3,4) +rho(3,3) +rho(3,2) +I*rho(3,1)
                             -I*rho(2,4) +rho(2,3) +rho(2,2) +I*rho(2,1)
                             -rho(1,4) -I*rho(1,3) -I*rho(1,2) +rho(1,1));

        rl[k] = 0.25 * creal( rho(4,4) -I*rho(4,3) +I*rho(4,2) +rho(4,1)
                             +I*rho(3,4) +rho(3,3) -rho(3,2) +I*rho(3,1)
                             -I*rho(2,4) -rho(2,3) +rho(2,2) -I*rho(2,1)
                             +rho(1,4) -I*rho(1,3) +I*rho(1,2) +rho(1,1));

        lh[k] = 0.50 * creal(rho(4,4) -I*rho(4,2) +I*rho(2,4) +rho(2,2));
        lv[k] = 0.50 * creal(rho(3,3) -I*rho(3,1) +I*rho(1,3) +rho(1,1));

        ld[k] = 0.25 * creal( rho(4,4) +rho(4,3) -I*rho(4,2) -I*rho(4,1)
                             +rho(3,4) +rho(3,3) -I*rho(3,2) -I*rho(3,1)
                             +I*rho(2,4) +I*rho(2,3) +rho(2,2) +rho(2,1)
                             +I*rho(1,4) +I*rho(1,3) +rho(1,2) +rho(1,1));

        la[k] = 0.25 * creal( rho(4,4) -rho(4,3) -I*rho(4,2) +I*rho(4,1)
                             -rho(3,4) +rho(3,3) +I*rho(3,2) -I*rho(3,1)
                             +I*rho(2,4) -I*rho(2,3) +rho(2,2) -rho(2,1)
                             -I*rho(1,4) +I*rho(1,3) -rho(1,2) +rho(1,1));

        lr[k] = 0.25 * creal( rho(4,4) +I*rho(4,3) -I*rho(4,2) +rho(4,1)
                             -I*rho(3,4) +rho(3,3) -rho(3,2) -I*rho(3,1)
                             +I*rho(2,4) -rho(2,3) +rho(2,2) +I*rho(2,1)
                             +rho(1,4) +I*rho(1,3) -I*rho(1,2) +rho(1,1));

        ll[k] = 0.25 * creal( rho(4,4) -I*rho(4,3) -I*rho(4,2) -rho(4,1)
                             +I*rho(3,4) +rho(3,3) +rho(3,2) -I*rho(3,1)
                             +I*rho(2,4) +rho(2,3) +rho(2,2) -I*rho(2,1)
                             -rho(1,4) +I*rho(1,3) +I*rho(1,2) +rho(1,1));


        PTrho[ 0] = exp(GAM_R*t) * creal(rho(4,4));
        PTrho[ 1] = exp(GAM_R*t) * cimag(rho(4,4));
        PTrho[ 2] = exp(GAM_R*t) * creal(rho(3,4));
        PTrho[ 3] = exp(GAM_R*t) * cimag(rho(3,4));
        PTrho[ 4] = exp(GAM_R*t) * creal(rho(4,2));
        PTrho[ 5] = exp(GAM_R*t) * cimag(rho(4,2));
        PTrho[ 6] = exp(GAM_R*t) * creal(rho(3,2));
        PTrho[ 7] = exp(GAM_R*t) * cimag(rho(3,2));

        PTrho[ 8] = exp(GAM_R*t) * creal(rho(4,3));
        PTrho[ 9] = exp(GAM_R*t) * cimag(rho(4,3));
        PTrho[10] = exp(GAM_R*t) * creal(rho(3,3));
        PTrho[11] = exp(GAM_R*t) * cimag(rho(3,3));
        PTrho[12] = exp(GAM_R*t) * creal(rho(4,1));
        PTrho[13] = exp(GAM_R*t) * cimag(rho(4,1));
        PTrho[14] = exp(GAM_R*t) * creal(rho(3,1));
        PTrho[15] = exp(GAM_R*t) * cimag(rho(3,1));

        PTrho[16] = exp(GAM_R*t) * creal(rho(2,4));
        PTrho[17] = exp(GAM_R*t) * cimag(rho(2,4));
        PTrho[18] = exp(GAM_R*t) * creal(rho(1,4));
        PTrho[19] = exp(GAM_R*t) * cimag(rho(1,4));
        PTrho[20] = exp(GAM_R*t) * creal(rho(2,2));
        PTrho[21] = exp(GAM_R*t) * cimag(rho(2,2));
        PTrho[22] = exp(GAM_R*t) * creal(rho(1,2));
        PTrho[23] = exp(GAM_R*t) * cimag(rho(1,2));

        PTrho[24] = exp(GAM_R*t) * creal(rho(2,3));
        PTrho[25] = exp(GAM_R*t) * cimag(rho(2,3));
        PTrho[26] = exp(GAM_R*t) * creal(rho(1,3));
        PTrho[27] = exp(GAM_R*t) * cimag(rho(1,3));
        PTrho[28] = exp(GAM_R*t) * creal(rho(2,1));
        PTrho[29] = exp(GAM_R*t) * cimag(rho(2,1));
        PTrho[30] = exp(GAM_R*t) * creal(rho(1,1));
        PTrho[31] = exp(GAM_R*t) * cimag(rho(1,1));

        //switch to basis <HH|, <H,V|, <V,H|, <V,V| (before it was the other way around)
        //PT = partially transposed


        m = gsl_matrix_complex_view_array (PTrho, 4, 4);

        gsl_eigen_herm(&m.matrix,eval,w);

        sum=0.;
        for(z=0;z<4;z++){

            value[z] = 0.0;} // wird auf 0 gesetzt für neuen Durchlauf der k-Schleife

        for(z=0;z<4;z++){

            value[z] = gsl_vector_get (eval, z);

            NEG[k] += fabs(value[z]) -value[z]; //only add the negative eigenvalue

            //sum+= fabs(value[i]) -value[i];

        }

        NEG[k] = 0.5 * NEG[k]; // the negative eigenvalue was added twice


          for(int l = k;l<k+200 && l<MAX_STEPS;l++){
            avePTrho[ 0][l] += exp(GAM_R*t) * creal(rho(4,4)) * delta_t;
            avePTrho[ 1][l] += exp(GAM_R*t) * cimag(rho(4,4)) * delta_t;
            avePTrho[ 2][l] += exp(GAM_R*t) * creal(rho(3,4)) * delta_t;
            avePTrho[ 3][l] += exp(GAM_R*t) * cimag(rho(3,4)) * delta_t;
            avePTrho[ 4][l] += exp(GAM_R*t) * creal(rho(4,2)) * delta_t;
            avePTrho[ 5][l] += exp(GAM_R*t) * cimag(rho(4,2)) * delta_t;
            avePTrho[ 6][l] += exp(GAM_R*t) * creal(rho(3,2)) * delta_t;
            avePTrho[ 7][l] += exp(GAM_R*t) * cimag(rho(3,2)) * delta_t;

            avePTrho[ 8][l] += exp(GAM_R*t) * creal(rho(4,3)) * delta_t;
            avePTrho[ 9][l] += exp(GAM_R*t) * cimag(rho(4,3)) * delta_t;
            avePTrho[10][l] += exp(GAM_R*t) * creal(rho(3,3)) * delta_t;
            avePTrho[11][l] += exp(GAM_R*t) * cimag(rho(3,3)) * delta_t;
            avePTrho[12][l] += exp(GAM_R*t) * creal(rho(4,1)) * delta_t;
            avePTrho[13][l] += exp(GAM_R*t) * cimag(rho(4,1)) * delta_t;
            avePTrho[14][l] += exp(GAM_R*t) * creal(rho(3,1)) * delta_t;
            avePTrho[15][l] += exp(GAM_R*t) * cimag(rho(3,1)) * delta_t;

            avePTrho[16][l] += exp(GAM_R*t) * creal(rho(2,4)) * delta_t;
            avePTrho[17][l] += exp(GAM_R*t) * cimag(rho(2,4)) * delta_t;
            avePTrho[18][l] += exp(GAM_R*t) * creal(rho(1,4)) * delta_t;
            avePTrho[19][l] += exp(GAM_R*t) * cimag(rho(1,4)) * delta_t;
            avePTrho[20][l] += exp(GAM_R*t) * creal(rho(2,2)) * delta_t;
            avePTrho[21][l] += exp(GAM_R*t) * cimag(rho(2,2)) * delta_t;
            avePTrho[22][l] += exp(GAM_R*t) * creal(rho(1,2)) * delta_t;
            avePTrho[23][l] += exp(GAM_R*t) * cimag(rho(1,2)) * delta_t;

            avePTrho[24][l] += exp(GAM_R*t) * creal(rho(2,3)) * delta_t;
            avePTrho[25][l] += exp(GAM_R*t) * cimag(rho(2,3)) * delta_t;
            avePTrho[26][l] += exp(GAM_R*t) * creal(rho(1,3)) * delta_t;
            avePTrho[27][l] += exp(GAM_R*t) * cimag(rho(1,3)) * delta_t;
            avePTrho[28][l] += exp(GAM_R*t) * creal(rho(2,1)) * delta_t;
            avePTrho[29][l] += exp(GAM_R*t) * cimag(rho(2,1)) * delta_t;
            avePTrho[30][l] += exp(GAM_R*t) * creal(rho(1,1)) * delta_t;
            avePTrho[31][l] += exp(GAM_R*t) * cimag(rho(1,1)) * delta_t;
        }

        if(t>100){

            double col[32] = {0.0};
            for(int i=0;i<32;i++){
              col[i] = avePTrho[i][k];
            }
            m = gsl_matrix_complex_view_array (col, 4, 4);

            gsl_eigen_herm(&m.matrix,eval,w);

            for(z=0;z<4;z++){

                ave_value[z][k] = gsl_vector_get (eval, z);

                ave_NEG[k] += fabs(ave_value[z][k]) -ave_value[z][k];

                //sum+= fabs(value[i]) -value[i];

            }
            //if(k>300&&k<305){printf("%f\t", ave_NEG[k]);}

            ave_NEG[k] = 0.5 * ave_NEG[k]*0.01;
        }


    // Print time-resolved Probabilities (numerics and analytics) and Negativity

        fprintf(f_pro_HX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",t*0.001,hh[k],hv[k],hd[k],ha[k],hr[k],hl[k]);
        fprintf(f_pro_VX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",t*0.001,vh[k],vv[k],vd[k],va[k],vr[k],vl[k]);
        fprintf(f_pro_DX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",t*0.001,dh[k],dv[k],dd[k],da[k],dr[k],dl[k]);
        fprintf(f_pro_AX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",t*0.001,ah[k],av[k],ad[k],aa[k],ar[k],al[k]);
        fprintf(f_pro_RX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",t*0.001,rh[k],rv[k],rd[k],ra[k],rr[k],rl[k]);
        fprintf(f_pro_LX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",t*0.001,lh[k],lv[k],ld[k],la[k],lr[k],ll[k]);

        fprintf(f_pro_ana_HX,"%g \t %g \t %g \t %g \t %g \t %g \t %g \n",t*0.001,P_ana(0,0,0,0,FSS,t),P_ana(0,0,M_PI,0,FSS,t),P_ana(0,0,0.5*M_PI,0,FSS,t),P_ana(0,0,0.5*M_PI,M_PI,FSS,t),P_ana(0,0,0.5*M_PI,0.5*M_PI,FSS,t),P_ana(0,0,0.5*M_PI,1.5*M_PI,FSS,t));

        fprintf(f_pro_ana_VX,"%g \t %g \t %g \t %g \t %g \t %g \t %g \n",t*0.001,P_ana(M_PI,0,0,0,FSS,t),P_ana(M_PI,0,M_PI,0,FSS,t),P_ana(M_PI,0,0.5*M_PI,0,FSS,t),P_ana(M_PI,0,0.5*M_PI,M_PI,FSS,t),P_ana(M_PI,0,0.5*M_PI,0.5*M_PI,FSS,t),P_ana(M_PI,0,0.5*M_PI,1.5*M_PI,FSS,t));

        fprintf(f_pro_ana_DX,"%g \t %g \t %g \t %g \t %g \t %g \t %g \n",t*0.001,P_ana(0.5*M_PI,0,0,0,FSS,t),P_ana(0.5*M_PI,0,M_PI,0,FSS,t),P_ana(0.5*M_PI,0,0.5*M_PI,0,FSS,t),P_ana(0.5*M_PI,0,0.5*M_PI,M_PI,FSS,t),P_ana(0.5*M_PI,0,0.5*M_PI,0.5*M_PI,FSS,t),P_ana(0.5*M_PI,0,0.5*M_PI,1.5*M_PI,FSS,t));

        fprintf(f_pro_ana_AX,"%g \t %g \t %g \t %g \t %g \t %g \t %g \n",t*0.001,P_ana(0.5*M_PI,M_PI,0,0,FSS,t),P_ana(0.5*M_PI,M_PI,M_PI,0,FSS,t),P_ana(0.5*M_PI,M_PI,0.5*M_PI,0,FSS,t),P_ana(0.5*M_PI,M_PI,0.5*M_PI,M_PI,FSS,t),P_ana(0.5*M_PI,M_PI,0.5*M_PI,0.5*M_PI,FSS,t),P_ana(0.5*M_PI,M_PI,0.5*M_PI,1.5*M_PI,FSS,t));

        fprintf(f_pro_ana_RX,"%g \t %g \t %g \t %g \t %g \t %g \t %g \n",t*0.001,P_ana(0.5*M_PI,0.5*M_PI,0,0,FSS,t),P_ana(0.5*M_PI,0.5*M_PI,M_PI,0,FSS,t),P_ana(0.5*M_PI,0.5*M_PI,0.5*M_PI,0,FSS,t),P_ana(0.5*M_PI,0.5*M_PI,0.5*M_PI,M_PI,FSS,t),P_ana(0.5*M_PI,0.5*M_PI,0.5*M_PI,0.5*M_PI,FSS,t),P_ana(0.5*M_PI,0.5*M_PI,0.5*M_PI,1.5*M_PI,FSS,t));

        fprintf(f_pro_ana_LX,"%g \t %g \t %g \t %g \t %g \t %g \t %g \n",t*0.001,P_ana(0.5*M_PI,1.5*M_PI,0,0,FSS,t),P_ana(0.5*M_PI,1.5*M_PI,M_PI,0,FSS,t),P_ana(0.5*M_PI,1.5*M_PI,0.5*M_PI,0,FSS,t),P_ana(0.5*M_PI,1.5*M_PI,0.5*M_PI,M_PI,FSS,t),P_ana(0.5*M_PI,1.5*M_PI,0.5*M_PI,0.5*M_PI,FSS,t),P_ana(0.5*M_PI,1.5*M_PI,0.5*M_PI,1.5*M_PI,FSS,t));


        fprintf(f_neg,"%g \t %g \n",t*0.001,NEG[k]);

        if(t<=100){
          fprintf(f_av_neg, "%g \t %g \n",t*0.001,0.0);
        }
        if(t>100){
        fprintf(f_av_neg, "%g \t %g \n",t*0.001,ave_NEG[k]);
        }



    } // end of time integration

    gsl_eigen_herm_free(w);
    gsl_vector_free(eval);

    // Here the numerically calculated probabilities are convolved by the time resolution of the detector (Gaussian), RESOLUTION = FWHM



    double NORMALIZATION; double GAUSSIAN;double tau;double t_fold;

    double FOLDED_HH,FOLDED_HV,FOLDED_HD,FOLDED_HA,FOLDED_HR,FOLDED_HL;
    double FOLDED_VH,FOLDED_VV,FOLDED_VD,FOLDED_VA,FOLDED_VR,FOLDED_VL;
    double FOLDED_DH,FOLDED_DV,FOLDED_DD,FOLDED_DA,FOLDED_DR,FOLDED_DL;
    double FOLDED_AH,FOLDED_AV,FOLDED_AD,FOLDED_AA,FOLDED_AR,FOLDED_AL;
    double FOLDED_RH,FOLDED_RV,FOLDED_RD,FOLDED_RA,FOLDED_RR,FOLDED_RL;
    double FOLDED_LH,FOLDED_LV,FOLDED_LD,FOLDED_LA,FOLDED_LR,FOLDED_LL;

    //k_MAX = k_MAX+1+k;

    for(k=0;k<MAX_STEPS;k++)
    { if(k==0) {k++;}

        FOLDED_HH=0.;FOLDED_HV=0.;FOLDED_HD=0.;FOLDED_HA=0.;FOLDED_HR=0.;FOLDED_HL=0.;
        FOLDED_VH=0.;FOLDED_VV=0.;FOLDED_VD=0.;FOLDED_VA=0.;FOLDED_VR=0.;FOLDED_VL=0.;
        FOLDED_DH=0.;FOLDED_DV=0.;FOLDED_DD=0.;FOLDED_DA=0.;FOLDED_DR=0.;FOLDED_DL=0.;
        FOLDED_AH=0.;FOLDED_AV=0.;FOLDED_AD=0.;FOLDED_AA=0.;FOLDED_AR=0.;FOLDED_AL=0.;
        FOLDED_RH=0.;FOLDED_RV=0.;FOLDED_RD=0.;FOLDED_RA=0.;FOLDED_RR=0.;FOLDED_RL=0.;
        FOLDED_LH=0.;FOLDED_LV=0.;FOLDED_LD=0.;FOLDED_LA=0.;FOLDED_LR=0.;FOLDED_LL=0.;

        tau=k*delta_t;

        NORMALIZATION = 0.;

        for(i=0;i<MAX_STEPS;i++)
        { if(i==0) {i++;}

            t = delta_t * i;
            t_fold = tau-t;
            GAUSSIAN = exp(-(t_fold-0.0*RESOLUTION)*(t_fold-0.0*RESOLUTION)/(2.*RESOLUTION*RESOLUTION));
            NORMALIZATION += delta_t * GAUSSIAN;

            FOLDED_HH += delta_t * hh[i] * GAUSSIAN;
            FOLDED_HV += delta_t * hv[i] * GAUSSIAN;
            FOLDED_HD += delta_t * hd[i] * GAUSSIAN;
            FOLDED_HA += delta_t * ha[i] * GAUSSIAN;
            FOLDED_HR += delta_t * hr[i] * GAUSSIAN;
            FOLDED_HL += delta_t * hl[i] * GAUSSIAN;

            FOLDED_VH += delta_t * vh[i] * GAUSSIAN;
            FOLDED_VV += delta_t * vv[i] * GAUSSIAN;
            FOLDED_VD += delta_t * vd[i] * GAUSSIAN;
            FOLDED_VA += delta_t * va[i] * GAUSSIAN;
            FOLDED_VR += delta_t * vr[i] * GAUSSIAN;
            FOLDED_VL += delta_t * vl[i] * GAUSSIAN;

            FOLDED_DH += delta_t * dh[i] * GAUSSIAN;
            FOLDED_DV += delta_t * dv[i] * GAUSSIAN;
            FOLDED_DD += delta_t * dd[i] * GAUSSIAN;
            FOLDED_DA += delta_t * da[i] * GAUSSIAN;
            FOLDED_DR += delta_t * dr[i] * GAUSSIAN;
            FOLDED_DL += delta_t * dl[i] * GAUSSIAN;

            FOLDED_AH += delta_t * ah[i] * GAUSSIAN;
            FOLDED_AV += delta_t * av[i] * GAUSSIAN;
            FOLDED_AD += delta_t * ad[i] * GAUSSIAN;
            FOLDED_AA += delta_t * aa[i] * GAUSSIAN;
            FOLDED_AR += delta_t * ar[i] * GAUSSIAN;
            FOLDED_AL += delta_t * al[i] * GAUSSIAN;

            FOLDED_RH += delta_t * rh[i] * GAUSSIAN;
            FOLDED_RV += delta_t * rv[i] * GAUSSIAN;
            FOLDED_RD += delta_t * rd[i] * GAUSSIAN;
            FOLDED_RA += delta_t * ra[i] * GAUSSIAN;
            FOLDED_RR += delta_t * rr[i] * GAUSSIAN;
            FOLDED_RL += delta_t * rl[i] * GAUSSIAN;

            FOLDED_LH += delta_t * lh[i] * GAUSSIAN;
            FOLDED_LV += delta_t * lv[i] * GAUSSIAN;
            FOLDED_LD += delta_t * ld[i] * GAUSSIAN;
            FOLDED_LA += delta_t * la[i] * GAUSSIAN;
            FOLDED_LR += delta_t * lr[i] * GAUSSIAN;
            FOLDED_LL += delta_t * ll[i] * GAUSSIAN;


            NORMALIZATION += delta_t * GAUSSIAN;
        }

        FOLDED_HH /= NORMALIZATION;
        FOLDED_HV /= NORMALIZATION;
        FOLDED_HD /= NORMALIZATION;
        FOLDED_HA /= NORMALIZATION;
        FOLDED_HR /= NORMALIZATION;
        FOLDED_HL /= NORMALIZATION;

        FOLDED_VH  /= NORMALIZATION;
        FOLDED_VV  /= NORMALIZATION;
        FOLDED_VD  /= NORMALIZATION;
        FOLDED_VA  /= NORMALIZATION;
        FOLDED_VR  /= NORMALIZATION;
        FOLDED_VL  /= NORMALIZATION;

        FOLDED_DH  /= NORMALIZATION;
        FOLDED_DV  /= NORMALIZATION;
        FOLDED_DD  /= NORMALIZATION;
        FOLDED_DA  /= NORMALIZATION;
        FOLDED_DR  /= NORMALIZATION;
        FOLDED_DL  /= NORMALIZATION;

        FOLDED_AH  /= NORMALIZATION;
        FOLDED_AV  /= NORMALIZATION;
        FOLDED_AD  /= NORMALIZATION;
        FOLDED_AA  /= NORMALIZATION;
        FOLDED_AR  /= NORMALIZATION;
        FOLDED_AL  /= NORMALIZATION;

        FOLDED_RH  /= NORMALIZATION;
        FOLDED_RV  /= NORMALIZATION;
        FOLDED_RD  /= NORMALIZATION;
        FOLDED_RA  /= NORMALIZATION;
        FOLDED_RR  /= NORMALIZATION;
        FOLDED_RL  /= NORMALIZATION;

        FOLDED_LH  /= NORMALIZATION;
        FOLDED_LV  /= NORMALIZATION;
        FOLDED_LD  /= NORMALIZATION;
        FOLDED_LA  /= NORMALIZATION;
        FOLDED_LR  /= NORMALIZATION;
        FOLDED_LL  /= NORMALIZATION;



        fprintf(f_pro_fold_HX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",tau*0.001,FOLDED_HH,FOLDED_HV,FOLDED_HD,FOLDED_HA,FOLDED_HR,FOLDED_HL);
        fprintf(f_pro_fold_VX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",tau*0.001,FOLDED_VH,FOLDED_VV,FOLDED_VD,FOLDED_VA,FOLDED_VR,FOLDED_VL);
        fprintf(f_pro_fold_DX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",tau*0.001,FOLDED_DH,FOLDED_DV,FOLDED_DD,FOLDED_DA,FOLDED_DR,FOLDED_DL);
        fprintf(f_pro_fold_AX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",tau*0.001,FOLDED_AH,FOLDED_AV,FOLDED_AD,FOLDED_AA,FOLDED_AR,FOLDED_AL);
        fprintf(f_pro_fold_RX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",tau*0.001,FOLDED_RH,FOLDED_RV,FOLDED_RD,FOLDED_RA,FOLDED_RR,FOLDED_RL);
        fprintf(f_pro_fold_LX,"%6g \t %8.6g \t %8.6g \t%8.6g \t%8.6g \t%8.6g \t%8.6g \n",tau*0.001,FOLDED_LH,FOLDED_LV,FOLDED_LD,FOLDED_LA,FOLDED_LR,FOLDED_LL);
        
        fprintf(fit_HH, "%6g \t %8.6g \n", tau*0.001, FOLDED_HH);
        fprintf(fit_HV, "%6g \t %8.6g \n", tau*0.001, FOLDED_HV);
        fprintf(fit_HD, "%6g \t %8.6g \n", tau*0.001, FOLDED_HD);
        fprintf(fit_HR, "%6g \t %8.6g \n", tau*0.001, FOLDED_HR);
        fprintf(fit_VH, "%6g \t %8.6g \n", tau*0.001, FOLDED_VH);
        fprintf(fit_VV, "%6g \t %8.6g \n", tau*0.001, FOLDED_VV);
        fprintf(fit_VD, "%6g \t %8.6g \n", tau*0.001, FOLDED_VD);
        fprintf(fit_VR, "%6g \t %8.6g \n", tau*0.001, FOLDED_VR);
        fprintf(fit_DH, "%6g \t %8.6g \n", tau*0.001, FOLDED_DH);
        fprintf(fit_DV, "%6g \t %8.6g \n", tau*0.001, FOLDED_DV);
        fprintf(fit_DD, "%6g \t %8.6g \n", tau*0.001, FOLDED_DD);
        fprintf(fit_DR, "%6g \t %8.6g \n", tau*0.001, FOLDED_DR);
        fprintf(fit_RH, "%6g \t %8.6g \n", tau*0.001, FOLDED_RH);
        fprintf(fit_RV, "%6g \t %8.6g \n", tau*0.001, FOLDED_RV);
        fprintf(fit_RD, "%6g \t %8.6g \n", tau*0.001, FOLDED_RD);
        fprintf(fit_RR, "%6g \t %8.6g \n", tau*0.001, FOLDED_RR);


    } // end of convolution








    // free given arrays
    free(temp);free(derivates); /* free(derivates_stst); */free(k1);free(k2);free(OUTPUT); //free(correlation);

    // close files

    fclose(f_pro_HX);fclose(f_pro_VX);fclose(f_pro_DX);fclose(f_pro_AX);fclose(f_pro_RX);fclose(f_pro_LX);
    fclose(f_pro_fold_HX);fclose(f_pro_fold_VX);fclose(f_pro_fold_DX);fclose(f_pro_fold_AX);fclose(f_pro_fold_RX);fclose(f_pro_fold_LX);

    fclose(f_neg);fclose(f_av_neg);
    
    fclose(fit_HH);fclose(fit_HV);fclose(fit_HD);fclose(fit_HR);fclose(fit_VH);fclose(fit_VV);fclose(fit_VD);fclose(fit_VR);fclose(fit_DH);fclose(fit_DV);fclose(fit_DD);fclose(fit_DR);fclose(fit_RH);fclose(fit_RV);fclose(fit_RD);fclose(fit_RR);


    return 0;
}
