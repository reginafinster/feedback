#include "itensor/all.h"

using namespace itensor;

double TLS_occupation(ITensor A)
{
    //search for Index in Tensor (should be system)
    Index s=findIndex(A, "sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(2),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real(); 
}

double TLS_norm(ITensor A)
{
    return eltC( dag(A)*A).real(); 
}


void SWAP_FORWARD(MPS& psi,const std::vector<Index>& bin, int to, int from, double cutoff)
{
    ITensor SWAP,U,S,V;
    Index iFB;
    for(int k=from;k<to;k++) // swap feedback bin next to tls bin
    {
            SWAP = psi.A(k)*psi.A(k+1); 
            iFB = findIndex(psi.A(k+1),"bath");
            if ( order(psi(k))==2 ) U=ITensor(iFB);
                              else  U=ITensor(iFB,commonIndex(psi.A(k),psi.A(k-1)));  
            svd(SWAP,U,S,V,{"Cutoff=",cutoff});
            psi.setA(k,U); 
            psi.setA(k+1,V*S); // orthoCenter wanders to system bin
    }
} 

void SWAP_BACKWARD(MPS& psi,const std::vector<Index>& bin, int from, int to, double cutoff)
{
    ITensor SWAP,U,S,V;
    Index iFB;
     for(int k=from;k>to;k--)
     {      
            SWAP = psi.A(k)*psi.A(k-1);
            iFB = findIndex(psi.A(k),"bath"); 
            if ( order(psi(k-1)) == 2 )  U=ITensor(iFB); // order(psi(k-1))==2, no commonIndex
                                  else   U=ITensor(iFB,commonIndex(psi.A(k-1),psi.A(k-2)));
            svd(SWAP,U,S,V,{"Cutoff=",cutoff}); 
            // -- lasse Orthocenter beim vorletzten Bin, der im naechsten Schritt wieder zum Systembin geswappt wird
            if (k-1==to)   { psi.setA(k-1,U);       psi.setA(k,V*S); }
            else           { psi.setA(k-1,U*S);     psi.setA(k,V);   }
     }
} 

//void MANIPULATE(Complex test[], double real, double imag){ test[0]=real+Cplx_i*imag;  }


int main(int argc, char* argv[])
    { 
        
time_t curtime; struct tm *loctime; curtime = time(NULL); loctime = localtime (&curtime); 
int hour = loctime -> tm_hour;int minute = loctime -> tm_min;int second = loctime -> tm_sec;
int year = 1900+loctime -> tm_year;int month = loctime -> tm_mon;int day  = loctime -> tm_mday;
printf("Date: %d.%d.%.d -- Time: %d:%d:%d  \n",day,month+1,year,hour,minute,second);  
        
    if(argc != 2) 
      { 
      //reminds us to give an input file if we forget
      printfln("Usage: %s inputfile",argv[0]); 
      return 0; 
      }
      
    auto input = InputGroup(argv[1],"input");
    //timestep
    Real dt = input.getReal("time_step");
    //end of integration
    int t_end = input.getInt("time_end");
    //number of steps
    //coherent pumping strengths
    Real Omega = input.getReal("Omega_TLS",0.);
    //TLS decay rate
    Real Gamma = input.getReal("Gamma",0.);
    //dimension of local Hilbert space of each bin
    int Nbin = input.getInt("Nbin",4); 
    //cutoff of schmidt values
    Real cutoff = input.getReal("svdcutoff");
    //maximal number of schmidtvalues 
    int maxm = input.getInt("maxnumberofSV");    
    int Nfb = input.getInt("feedback_time");
    int ttotal = t_end + Nfb;
    int SHOW_PROGRESS =1;
    int SAVE_EVERY_STEP = 1; 
    if (ttotal>1000) { SHOW_PROGRESS=100; SAVE_EVERY_STEP=ttotal/10; }
    
    Real init_22 = input.getReal("init_excited_state"); 
    Real init_11 = input.getReal("init_ground_state"); 
    Real phi = input.getReal("phi"); 
    Real fb_on = input.getReal("fb_on");
    
    int SpectrumSteps = input.getInt("SpectrumSteps");
    Real SpectrumIntervall = input.getReal("SpectrumIntervall");    
    
    printf("TLS initial state: (22)=%.2f -- (11)=%.2f \n",init_22,init_11);
    
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- SETUP THE MPS ----------------------------------------------
    // ----------------------------------------------------------------------------------
    // create the physical indices
    
    auto bin = std::vector<Index>((int)ttotal+1);
    int tls=0; bin.at(tls) = Index(2,"sys"); // bin(0) sys throughout the code
    for(int j = 1; j <= ttotal; ++j){ bin.at(j) = Index(Nbin,"bath"); } // ttotal steps

    auto binlink = std::vector<Index>((int)ttotal+1);
    for(int j = 1; j <= ttotal; ++j) { binlink.at(j) = Index(1,"Link");   }
    MPS psi=MPS(ttotal+1);
    ITensor MPSTensor = ITensor();

    int now  = Nfb+1;
    int past = now-Nfb;
    
    MPSTensor = ITensor(bin[past],binlink[1]);
    MPSTensor.set(bin[past](1),binlink[1](1),1.); 
    psi.setA((past),MPSTensor);

    for(int j = 2; j<=Nfb; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j-1],binlink[j]);
     MPSTensor.set(bin[j](1),binlink[j-1](1),binlink[j](1),1.); // Level besetzt im TLS
     psi.setA((j),MPSTensor);
    }   
    // until the system bin psi(i) contain bin[i]
    // at MPS slot Nfb the tls tensor is put
    MPSTensor = ITensor(bin[tls],binlink[Nfb],binlink[Nfb+1]);
    MPSTensor.set(bin[tls](2),binlink[Nfb](1),binlink[Nfb+1](1),init_22); // Level besetzt im TLS
    MPSTensor.set(bin[tls](1),binlink[Nfb](1),binlink[Nfb+1](1),init_11); // Level besetzt im TLS
    psi.setA((now),MPSTensor);

    for(int j = Nfb+1; j<ttotal; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j],binlink[j+1]);
     MPSTensor.set(bin[j](1),binlink[j](1),binlink[j+1](1),1.); // Level besetzt im TLS
     psi.setA((j+1),MPSTensor);
    }
    // the last slot of MPS is filled with a tensor with only links to the left
    MPSTensor = ITensor(bin[ttotal],binlink[ttotal]);
    MPSTensor.set(bin[ttotal](1),binlink[ttotal](1),1.); // Level besetzt im TLS
    psi.setA((ttotal+1),MPSTensor);
    
    //for(int l=0;l<=ttotal+1;l++) {printf("%i:",l); PrintData(psi(l));} printf("Done!\n");
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------END MPS SETUP -----------------------------------------------
    // ----------------------------------------------------------------------------------

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------U_EVO SETUP -------------------------------------------------
    // ----------------------------------------------------------------------------------
    auto H_sys =ITensor(bin[tls],prime(bin[tls]));
    // interaction with futurebin
    auto H_dis = ITensor(bin[tls],prime(bin[tls]),bin[now],prime(bin[now]));
    // interaction with the pastbin
    auto H_fb    = ITensor(bin[tls],prime(bin[tls]),bin[past],prime(bin[past]));
        
    //pumped TLS 
    H_sys.set(bin[tls](1),prime(bin[tls](2)),-Cplx_i*dt*Omega);
    H_sys.set(bin[tls](2),prime(bin[tls](1)),-Cplx_i*dt*Omega);
    
    //H_int with photon creation/destruction in bin k
   for(int phot=1;phot<Nbin;phot++)
   {            
    H_dis.set(bin[tls](2),prime(bin[tls](1)),bin[now](phot  ) ,prime(bin[now](phot+1)) ,( 1.)*sqrt(dt)*Gamma*sqrt(phot));
    H_dis.set(bin[tls](1),prime(bin[tls](2)),bin[now](phot+1) ,prime(bin[now](phot  )) ,(-1.)*sqrt(dt)*Gamma*sqrt(phot));        
    H_fb.set(bin[tls](2) ,prime(bin[tls](1)),bin[past](phot  ),prime(bin[past](phot+1)),(-1.)*sqrt(dt)*Gamma*exp(-Cplx_i*phi)*sqrt(phot));
    H_fb.set(bin[tls](1) ,prime(bin[tls](2)),bin[past](phot+1),prime(bin[past](phot  )),( 1.)*sqrt(dt)*Gamma*exp( Cplx_i*phi)*sqrt(phot));
   }//phot

    auto H_int =  H_dis * delta(bin[past],prime(bin[past])) 
                + H_fb  * delta(bin[now],prime(bin[now]))    ;
    
    auto H_fb_1  = delta(bin[now],prime(bin[now])) * H_sys * delta(bin[past],prime(bin[past])) + H_int; 
    auto H_fb_2  = (1./2.)  * mapPrime(H_fb_1*prime(H_fb_1),2,1);
    auto H_fb_3  = (1./3.)  * mapPrime(H_fb_1*prime(H_fb_2),2,1);
    auto H_fb_4  = (1./4.)  * mapPrime(H_fb_1*prime(H_fb_3),2,1);
    auto H_fb_5  = (1./5.)  * mapPrime(H_fb_1*prime(H_fb_4),2,1);
    auto H_fb_6  = (1./6.)  * mapPrime(H_fb_1*prime(H_fb_5),2,1);
    auto H_fb_7  = (1./7.)  * mapPrime(H_fb_1*prime(H_fb_6),2,1);
    auto H_fb_8  = (1./8.)  * mapPrime(H_fb_1*prime(H_fb_7),2,1);
    auto H_fb_9  = (1./9.)  * mapPrime(H_fb_1*prime(H_fb_8),2,1);
    auto H_fb_10 = (1./10.) * mapPrime(H_fb_1*prime(H_fb_9),2,1);

    auto delta_temp = ITensor(bin[tls],prime(bin[tls]));
         delta_temp.set(bin[tls](1),prime(bin[tls](1)),1.);
         delta_temp.set(bin[tls](2),prime(bin[tls](2)),1.);
    
    auto U_evo =   delta(bin[now],prime(bin[now])) * delta_temp * delta(bin[past],prime(bin[past])) 
                 + H_fb_1 + H_fb_2 + H_fb_3 + H_fb_4 + H_fb_5 + H_fb_6 + H_fb_7  + H_fb_8 + H_fb_9 + H_fb_10 ;
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- END U_EVO SETUP --------------------------------------------
    // ----------------------------------------------------------------------------------

    
    ITensor U,S,V,SWAP,temp;
    int j=1;
    
    double cv_norm,cv_pop;
    
    // das bild im MPS sieht aus [pastbin][feedbackbin][systembin][currentbin]
    Index iFB,iSB,iCB,iPB; // index fuer feedback, system, und current bin 
    Index iFBPB; // index fuer feedback bin und past bin
    // zum zeit berechnen
    FILE * file;
    char FILE_NAME[2048+1024];
    snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f_phi_%.2f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.,phi);
    file = fopen(FILE_NAME,"w");
    fprintf(file,"#time \tS+S- \t norm \n"); fprintf(file,"# %02d_%02d_%02d_%02d_%02d_%02d \n",year,month+1,day,hour,minute,second);
    
    int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
    int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
    double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

for(int i=now;i<=ttotal;i++)            
{   
     if ( i % SAVE_EVERY_STEP == 0) {fprintf(file,"%.10f \t %.10f \t %10.f \n",(i-Nfb-1)*dt,cv_pop,cv_norm); fflush(file);} 
     // -------------------- SWAP bin in psi(i-Nfb) to psi(i-1) --------------------------------------------
     SWAP_FORWARD(psi,bin,now-1,past,cutoff); // moves ortoCenter in now-1 // does nothing if now-1-past=1
     // -----------------------------------------------------------------------------------------------
     // define the links for feedback bin, sys bin and current bin
     if (i-2>0) iFBPB = commonIndex(psi.A(i-1),psi.A(i-2)); // exclude case Nfb=1 for i=1
     iFB = findIndex(psi.A(i-1),"bath"); 
     iSB = findIndex(psi.A(i),"sys");
     iCB = findIndex(psi.A(i+1),"bath");
     // now apply the mpo
     temp=noPrime(U_evo * psi.A(i-1) * psi.A(i) * psi.A(i+1)); 
     U=ITensor(iCB,iFB,iFBPB); 
     svd(temp,U,S,V,{"Cutoff=",cutoff});
     // tls in the next mps slot 
     psi.setA(i+1,V); 
     cv_norm = TLS_norm(V*S); 
     cv_pop  = TLS_occupation(V*S); 
     // now factorize feedback bin and current bin
     temp=U*S; 
     U=ITensor(iFB,iFBPB); 
     svd(temp,U,S,V,{"Cutoff=",cutoff});
     if (Nfb>1) { psi.setA(i,V  ); psi.setA(i-1,U*S);} // if feedback length 1 no swap is necessary
     else       { psi.setA(i,V*S); psi.setA(i-1,U  );} // so the orthoCenter stays next to system
     // -------------------- SWAP bin in psi(i-1) to psi(i-Nfb) ---------------------------------------------
     SWAP_BACKWARD(psi,bin,now-1,past,cutoff); // keeps ortoCenter in past+1 // does nothing if now-1-past=1
     // -----------------------------------------------------------------------------------------------
     if (i<ttotal)
     {    
     // define new evolution operator 
     U_evo = delta(bin[now],bin[now+1])   * U_evo * delta(prime(bin[now] ),prime(bin[now+1])) ;
     U_evo = delta(bin[past],bin[past+1]) * U_evo * delta(prime(bin[past]),prime(bin[past+1]))   ;
     now++;past++;
     }
     // show the calculation status every percent
     if ( (i % SHOW_PROGRESS ) == 0)  
       { 
         curtime = time(NULL); loctime = localtime (&curtime); 
         ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
         time_now = esecond + 60*eminute+60*60*ehour;  
         printf("Step %i of %i -- pop=%.10f -- norm=%.10f -- Time elapsed=%.2f sec. \n",i,ttotal,cv_pop,cv_norm,time_now-time_before);
         time_before = time_now;
       }     
}//end timeloop
printf("End of loop\n");
// give information how long it took
curtime = time(NULL); loctime = localtime (&curtime); 
ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
eyear = 1900+loctime -> tm_year; emonth = loctime -> tm_mon; eday  = loctime -> tm_mday;
printf("End of calculation: %d.%d.%.d -- Time: %d:%d:%d  \n",eday,emonth+1,eyear,ehour,eminute,esecond);  
fprintf(file,"# end of calculation %02d_%02d_%02d_%02d_%02d_%02d \n",eyear,emonth+1,eday,ehour,eminute,esecond);
printf("Time elapsed: %d.%d.%.d -- Time: %dh:%dmin:%dsec \n",eday-day,emonth-month,eyear-year,ehour-hour,eminute-minute,esecond-second);  
fprintf(file,"# time difference %02d_%02d_%02d_%02d_%02d_%02d \n",eyear-year,emonth-month,eday-day,ehour-hour,eminute-minute,esecond-second);
// close the output file pointer    
fclose(file);

// ----------------------------------------------  time evolution end ------------------------------
// calculate observables

// PrintData(psi(ttotal+1)); // here is the system bin
//PrintData(psi(ttotal-Nfb+1)); // this is the next feedback bin with orthoCenter, therefore normalized
//temp = dag(psi(ttotal-Nfb+1))*psi(ttotal-Nfb+1); // this is the next feedback bin with OrthoCenter 
//PrintData(temp); 

//swap orthoCenter to the previous one
temp = psi(ttotal-Nfb+1)*psi(ttotal-Nfb);
U = ITensor(bin[ttotal-Nfb+1],commonIndex(psi(ttotal-Nfb+1),psi(ttotal-Nfb+2)));
svd(temp,U,S,V);
psi.setA(ttotal-Nfb+1,U);
psi.setA(ttotal-Nfb,V*S);
// check whether psi(ttotal-Nfb) is now normalized
temp = dag(psi(ttotal-Nfb))*psi(ttotal-Nfb);
printf("Check if psi with reference bin is normalized:\n"); PrintData(temp);
 // now the first not interacting with the system anymore is normalized
int RefBin=ttotal-Nfb;

snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_BathCorrelations_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f_phi_%.2f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.,phi);
file = fopen(FILE_NAME,"w");

// now we have a steady state and can calculate the spectrum <b^\dg(t)b(t-\tau)>
// setup the bath correlation array
auto bath_g1 = std::vector<Complex>(SpectrumSteps+1);
auto bath_g2= std::vector<Complex>(SpectrumSteps+1);
auto bath_pop = std::vector<Complex>(SpectrumSteps+1);
ITensor Tg1;
ITensor Tg2;
ITensor Tb_pop;
    
// now define the flip operator for RefBin
Index s=findIndex(psi(RefBin),"bath"); //Print(s);
ITensor Sp = ITensor(s,prime(s));
// watch of for the maximum photon number in the tensor
   for(int phot=1;phot<Nbin;phot++) Sp.set(s(phot+1),prime(s(phot)),sqrt(phot)*1.); //PrintData(Sp);
    
// the first entry is the self-correlation of the reference bin
Tg1= dag( noPrime(psi(RefBin)*Sp) )*noPrime(Sp*psi(RefBin));
Tg2= dag( noPrime(Sp*noPrime(Sp*psi(RefBin))) ) * noPrime(Sp*noPrime(Sp*psi(RefBin)));
Tb_pop = dag( noPrime(Sp*psi(RefBin)) ) * noPrime(Sp*psi(RefBin));
bath_g1[0]=eltC(Tg1);
bath_g2[0]=eltC(Tg2);
bath_pop[0]=eltC(Tb_pop);
int counter=0; // because easier to count
printf("bath_g1[%i]-|bc_inf|^2=%.10f+I%.10f -- g2=%.10f -- bath_pop[0]=%.10f  \n",counter,bath_g1[counter].real(),bath_g1[counter].imag(),bath_g2[counter].real(),bath_pop[counter].real());
fprintf(file,"%.10f \t %.10f \t %10f \t %.10f \n",dt*counter,bath_g1[counter].real(),bath_g1[counter].imag(),bath_g2[counter].real());

//Index for CorrelatedBin
Index t;
ITensor Sm;

for(int i=RefBin;i>RefBin-SpectrumSteps;i--)
{   
    t=findIndex(psi(i-1),"bath"); //Print(t);
    Sm = ITensor(t,prime(t));
    for(int phot=1;phot<Nbin;phot++) Sm.set(t(phot+1),prime(t(phot)),sqrt(phot)*1.);
    temp = psi(i)*psi(i-1);
    Tg1= dag( noPrime(temp*Sp) )*noPrime(Sm*temp); // (Sp*temp)^dg*(Sm*temp)=(B_R|psi>)^+ B_b|psi>=<B^+_R B_b>
    Tg2= dag( noPrime(Sm*noPrime(Sp*temp) )) * noPrime(Sm*noPrime(Sp*temp));
    Tb_pop = dag( noPrime(Sp*temp) ) * noPrime(Sp*temp);
    counter++;
    bath_g1[counter]=eltC(Tg1);
    bath_g2[counter]=eltC(Tg2);
    bath_pop[counter]=eltC(Tb_pop);
    printf("bath_g1[%i]-|bc_inf|^2=%.10f+I%.10f -- g2=%.10f -- (bath_pop)^2=%.10f\n",counter,bath_g1[counter].real(),bath_g1[counter].imag(),bath_g2[counter].real(),bath_pop[counter].real()*bath_pop[counter].real());
    fprintf(file,"%.10f \t %.10f \t %10f \t %.10f \n",dt*counter,bath_g1[counter].real(),bath_g1[counter].imag(),bath_g2[counter].real());

    U = ITensor(bin[RefBin],commonIndex(psi(i-1),psi(i-2))); //PrintData(U);
    svd(temp,U,S,Tg1);
    psi.setA(i-1,U*S);
    psi.setA(i,Tg1); 
}

Complex b_offset = bath_g1[counter]; // needed to calculate the spectrum without coherent part
double g2_norm = bath_g2[counter].real(); // needed for normalization // correct would be to use (bath_pop)^2, so check
fclose(file);    

// ---------------- output spectrum --------------------------------------
    
snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_Spectrum_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f_phi_%.2f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.,phi);
file = fopen(FILE_NAME,"w");

auto spectrum = std::vector<Complex>(SpectrumSteps+1);
double dw = SpectrumIntervall/SpectrumSteps; 
double om;
   
for(int i=0;i<=SpectrumSteps;i++)
{
     om = -0.5*(SpectrumIntervall)+i*dw;
     spectrum[i] = 0.;
     for(int b=0;b<=SpectrumSteps;b++) spectrum[i] += exp(-Cplx_i*om*b*dt)*(bath_g1[b]-b_offset);  
     fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \n",om,spectrum[i].real(),spectrum[i].imag(),bath_g2[i].real()/g2_norm);
}
fclose(file);

return 0;
}
