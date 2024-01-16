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
    Real t_end = input.getReal("time_end");
    //number of steps
    Real ttotal = t_end/dt;
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
    Real init_22 = input.getReal("init_excited_state"); 
    Real init_11 = input.getReal("init_ground_state"); 
    Real phi = input.getReal("phi"); 
    printf("TLS initial state: (22)=%.2f -- (11)=%.2f \n",init_22,init_11);
    if (Nfb == 0) Nfb = ttotal;
    //dimension of system bin
    int d = 2; 

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- SETUP THE MPS ----------------------------------------------
    // ----------------------------------------------------------------------------------
    auto bin = std::vector<Index>((int)ttotal+3);
    for(int j = 1; j <= ttotal+1; ++j){ bin.at(j) = Index(Nbin,"bath"); }
    int tls=ttotal+2;
    bin.at(tls) = Index(Nbin,"sys"); // bin(1) sys throughout the code

    auto binlink = std::vector<Index>((int)ttotal+1);
    for(int j = 1; j <= ttotal; ++j) { binlink.at(j) = Index(1,"Link");   }
    MPS psi=MPS(ttotal+1);
    ITensor MPSTensor = ITensor();

    // set up the mps 
    // first MPS slot is filled with a bin which has no link to the left
     MPSTensor = ITensor(bin[1],binlink[1]);
     MPSTensor.set(bin[1](1),binlink[1](1),1.); // Level besetzt im TLS
     psi.setA((1),MPSTensor);

    for(int j = 2; j<Nfb; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j-1],binlink[j]);
     MPSTensor.set(bin[j](1),binlink[j-1](1),binlink[j](1),1.); // Level besetzt im TLS
     psi.setA((j),MPSTensor);
    }   
    // at MPS slot Nfb the tls tensor is put
    MPSTensor = ITensor(bin[tls],binlink[Nfb-1],binlink[Nfb]);
    MPSTensor.set(bin[tls](2),binlink[Nfb-1](1),binlink[Nfb](1),init_22); // Level besetzt im TLS
    MPSTensor.set(bin[tls](1),binlink[Nfb-1](1),binlink[Nfb](1),init_11); // Level besetzt im TLS
    psi.setA((Nfb),MPSTensor);
    
    for(int j = Nfb; j<=ttotal-1; j++)
    {
     MPSTensor = ITensor(bin[j],binlink[j],binlink[j+1]);
     MPSTensor.set(bin[j](1),binlink[j](1),binlink[j+1](1),1.); // Level besetzt im TLS
     psi.setA((j+1),MPSTensor);
    }
    // the last slot of MPS is filled with a tensor with only links to the left
    MPSTensor = ITensor(bin[ttotal+1],binlink[ttotal]);
    MPSTensor.set(bin[ttotal+1](1),binlink[ttotal](1),1.); // Level besetzt im TLS
    psi.setA((ttotal+1),MPSTensor);
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------END MPS SETUP -----------------------------------------------
    // ----------------------------------------------------------------------------------

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // ---------------------U_EVO SETUP -------------------------------------------------
    // ----------------------------------------------------------------------------------
    auto U_sysfb =ITensor(bin[tls],prime(bin[tls]));
    // interaction with futurebin
    auto U_disfb = ITensor(bin[tls],prime(bin[tls]),bin[Nfb],prime(bin[Nfb]));
    // interaction with the pastbin
    auto U_fb    = ITensor(bin[tls],prime(bin[tls]),bin[1],prime(bin[1]));
        
    //pumped TLS 
    U_sysfb.set(bin[tls](1),prime(bin[tls](2)),-Cplx_i*dt*Omega);
    U_sysfb.set(bin[tls](2),prime(bin[tls](1)),-Cplx_i*dt*Omega);
    
    //H_int with photon creation/destruction in bin k
   for(int phot=1;phot<Nbin;phot++)
   {            
    U_disfb.set(bin[tls](2),prime(bin[tls](1)),bin[Nfb](phot  ),prime(bin[Nfb](phot+1)),( 1.)*sqrt(dt)*Gamma*sqrt(phot));
    U_disfb.set(bin[tls](1),prime(bin[tls](2)),bin[Nfb](phot+1),prime(bin[Nfb](phot  )),(-1.)*sqrt(dt)*Gamma*sqrt(phot));        
    U_fb.set(bin[tls](2),prime(bin[tls](1)),bin[1](phot  ),prime(bin[1](phot+1)),(-1.)*sqrt(dt)*Gamma*exp(-Cplx_i*phi)*sqrt(phot));
    U_fb.set(bin[tls](1),prime(bin[tls](2)),bin[1](phot+1),prime(bin[1](phot  )),( 1.)*sqrt(dt)*Gamma*exp(Cplx_i*phi)*sqrt(phot));
   }//phot

    auto U_int = U_disfb*delta(bin[1],prime(bin[1])) 
                +U_fb   *delta(bin[Nfb],prime(bin[Nfb]))    ;
    
    ITensor temp;            
    auto tempfb = ITensor(bin[tls],prime(bin[tls]));
    tempfb.set(bin[tls](1),prime(bin[tls](1)),1.);
    tempfb.set(bin[tls](2),prime(bin[tls](2)),1.);

    auto U_fb_1  = U_sysfb*delta(bin[Nfb],prime(bin[Nfb]))*delta(bin[1],prime(bin[1])) +U_int; 
    auto U_fb_2  = (1./2.)  * mapPrime(U_fb_1*prime(U_fb_1),2,1);
    auto U_fb_3  = (1./3.)  * mapPrime(U_fb_1*prime(U_fb_2),2,1);
    auto U_fb_4  = (1./4.)  * mapPrime(U_fb_1*prime(U_fb_3),2,1);
    auto U_fb_5  = (1./5.)  * mapPrime(U_fb_1*prime(U_fb_4),2,1);
    auto U_fb_6  = (1./6.)  * mapPrime(U_fb_1*prime(U_fb_5),2,1);
    auto U_fb_7  = (1./7.)  * mapPrime(U_fb_1*prime(U_fb_6),2,1);
    auto U_fb_8  = (1./8.)  * mapPrime(U_fb_1*prime(U_fb_7),2,1);
    auto U_fb_9  = (1./9.)  * mapPrime(U_fb_1*prime(U_fb_8),2,1);
    auto U_fb_10 = (1./10.) * mapPrime(U_fb_1*prime(U_fb_9),2,1);
    
    auto U_evofb =   delta(bin[Nfb],prime(bin[Nfb]))*tempfb*delta(bin[1],prime(bin[1])) + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7  + U_fb_8 + U_fb_9 + U_fb_10 ;
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // --------------------- END U_EVO SETUP --------------------------------------------
    // ----------------------------------------------------------------------------------

    
    ITensor U,S,V,W,SWAP;
    int j=1;
    
    double cv_norm,cv_pop;
    
// das bild im MPS sieht aus [pastbin][feedbackbin][systembin][currentbin]
Index iFB,iSB,iCB,iPB; // index fuer feedback, system, und current bin 
Index iFBPB; // index fuer feedback bin und past bin
// zum zeit berechnen
    FILE * file;
    char FILE_NAME[2048+1024];
    snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.);
    file = fopen(FILE_NAME,"w");
    fprintf(file,"#time \tS+S- \t norm \n");
    fprintf(file,"# %02d_%02d_%02d_%02d_%02d_%02d \n",year,month+1,day,hour,minute,second);
    fprintf(file,"%f \t %e \t %e \n",0.,init_22,0.0000);
    

int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

//for(int i=Nfb;i<ttotal;i++)            
for(int i=Nfb;i<ttotal;i++)            
{   //PrintData(psi(1));PrintData(psi(2));
    //PrintData(psi(i-3)); PrintData(psi(i-2)); PrintData(psi(i-1)); PrintData(psi(i)); PrintData(psi(i+1));
    // zuerst den feedbackbin nach i-1 schaufeln         
    for(int k=1;k<Nfb-1;k++) // immer eine Delayschlaufe zurueck
    {
            SWAP = psi.A(i-Nfb+k)*psi.A(i-Nfb+k+1); // if (k==Nfb-2) PrintData(SWAP);
            if (i-Nfb+k==1) U=ITensor(bin[i-Nfb+k+1]); // keine Linksverlinkung fuer psi.A(1)
            else U=ITensor(bin[i-Nfb+k+1],commonIndex(psi.A(i-Nfb+k),psi.A(i-Nfb+k-1))); // if (k==Nfb-2) PrintData(U); 
            svd(SWAP,U,S,V,{"Cutoff=",cutoff});//,{"LeftIndexName", binlink[i-Nfb+k].name(),"Cutoff",cutoff,"Maxm",maxm});
            psi.setA(i-Nfb+k,U); //if (k==Nfb-2) {PrintData(U); PrintData(V*S);}
            psi.setA(i-Nfb+k+1,S*V);
    }
      // ------ der Uebersicht wegen identifizierbare die zugehoerigen Links 
     iFB = findIndex(psi.A(i-1),"bath"); iFBPB = commonIndex(psi.A(i-1),psi.A(i-2));
     iSB = findIndex(psi.A(i),"sys");
     iCB = findIndex(psi.A(i+1),"bath");
     // ------ kontrahiere FB,SB, und CB und fuehre eine zweifache Schmidtzerlegung durch
 
     temp=noPrime(U_evofb * psi.A(i-1) * psi.A(i) * psi.A(i+1)); //PrintData(temp);
     U=ITensor(iCB,iFB,iFBPB); //PrintData(U);
     svd(temp,U,S,V,{"Cutoff=",cutoff});
     // Systembin nach zum naechsten Zeitschritt i+1
     psi.setA(i+1,V); //PrintData(V); PrintData(U);
     cv_norm = TLS_norm(V*S);
     cv_pop  = TLS_occupation(V*S);
     //printf("%.10f \t %.10f \t %10.f \n",i*dt,cv_pop,cv_norm);
     fprintf(file,"%.10f \t %.10f \t %10.f \n",(i-Nfb+1)*dt,cv_pop,cv_norm);
     //if ( (j % ((int)Nfb/10)) == 0) 
     //printf("time %.3f of %.3f -- pop=%.10f - norm=%.10f \n",j*dt,t_end*1.,cv_pop,cv_norm);
     // ---- Berechne den Erwartungswert und schiebe dafuer das Orthocenter zum Systembin
     //norm_psi=1.-TLS_norm(V*S,systembin);
     //fprintf(file,"%f\t%e\t%e\n",i*dt,TLS_occupation(V*S),norm(psi));
     // zerlege nun currentbin und feedbackbin und stelle feedbackbin nach links und gebe ihm das Orthocenter
     W=U*S; 
     U=ITensor(iFB,iFBPB); //PrintData(U);
     svd(W,U,S,V,{"Cutoff=",cutoff});
     // Feedbackbin nach i-1
     psi.setA(i-1,U*S); //PrintData(U);PrintData(V);
     // Currentbin nach i
     psi.setA(i,V);
     //psi.position(i);
     //PrintData(psi(i-2)); PrintData(psi(i-1)); PrintData(psi(i)); PrintData(psi(i+1)); 
     // ----- zurueck swappen, da fuer viele tau-intervalle das hinzuswappen zu lange dauern wuerde ---
     for(int k=1;k<Nfb-1;k++)
     {      
            SWAP = psi.A(i-1-k)*psi.A(i-k);
            if ( (i-k-1)==1 ) U=ITensor(bin[i-Nfb+1]); // kein zu vererbender Link
            else U=ITensor(bin[i-Nfb+1],commonIndex(psi.A(i-k-2),psi.A(i-k-1)));
            svd(SWAP,U,S,V,{"Cutoff=",cutoff});//,{"LeftIndexName", binlink[i-k-1].name(),"Cutoff",cutoff,"Maxm",maxm});
            // -- lasse Orthocenter beim vorletzten Bin, der im naechsten Schritt wieder zum Systembin geswappt wird
            if (k==Nfb-2) { psi.setA(i-1-k,U);       psi.setA(i-k,V*S); }
            else          { psi.setA(i-1-k,U*S);       psi.setA(i-k,V); }
     }
     U_evofb=U_evofb*delta(bin[i],bin[i+1])*delta(prime(bin[i]),prime(bin[i+1]))*delta(bin[i-Nfb+1],bin[i-Nfb+2])*delta(prime(bin[i-Nfb+1]),prime(bin[i-Nfb+2]));  
     // ----- gebe die Normabweichung aus und ueberpruefe ob MPS in kanonischer Form (maximal Tensor dritter Stufe)
     if ( (i % ((int)ttotal/100)) == 0) 
       { 
         curtime = time(NULL); loctime = localtime (&curtime); 
         ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
         time_now = esecond + 60*eminute+60*60*ehour;  
         printf("Feedback time %.3f of %.3f -- pop=%.10f -- norm=%.10f -- Time elapsed=%.2f sec. \n",i*dt,t_end*1.,cv_pop,cv_norm,time_now-time_before);
         time_before = time_now;
       }     
     fflush(file); // This will flush any pending fprintf output
    }//end timeloop

curtime = time(NULL); loctime = localtime (&curtime); 
ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
eyear = 1900+loctime -> tm_year; emonth = loctime -> tm_mon; eday  = loctime -> tm_mday;
printf("End of calculation: %d.%d.%.d -- Time: %d:%d:%d  \n",eday,emonth+1,eyear,ehour,eminute,esecond);  
fprintf(file,"# end of calculation %02d_%02d_%02d_%02d_%02d_%02d \n",eyear,emonth+1,eday,ehour,eminute,esecond);
printf("Time elapsed: %d.%d.%.d -- Time: %dh:%dmin:%dsec \n",eday-day,emonth-month,eyear-year,ehour-hour,eminute-minute,esecond-second);  
fprintf(file,"# time difference %02d_%02d_%02d_%02d_%02d_%02d \n",eyear-year,emonth-month,eday-day,ehour-hour,eminute-minute,esecond-second);
    
    fclose(file);
    return 0;
    }
