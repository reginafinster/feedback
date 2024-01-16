#include "itensor/all.h"

using namespace itensor;

double TLS_occupation(ITensor A, IndexType type)
{
    //search for Index in Tensor (should be system)
    Index s=findtype(A, type) ;       
    if(!s) printf("Error in calculating the TLS_occupation:\nBin has no Index of type %s!\n",type);
    ITensor Sm = ITensor(s,prime(s));
    Sm.set(s(2),prime(s(1)),1.);
    
    return ( dag(prime(A*Sm,s))*Sm*A ).real();
}
    
double TLS_norm(ITensor A, IndexType type)
{
    Index s=findtype(A, type) ;       
    if(!s) printf("Error in calculating the TLS_norm:\nBin has no Index of type %s!\n",type);
    
    return (dag(A) * A).real();
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
    int ttotal = input.getInt("ttotal",0); 
    //Real ttotal = (int)(t_end/dt);
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
    Real global_fak = input.getReal("fb_fak");
    printf("TLS initial state: (22)=%.2f -- (11)=%.2f \n",init_22,init_11);
    if (Nfb == 0) Nfb = ttotal;
    FILE * file;
    char FILE_NAME[2048+1024];
    snprintf(FILE_NAME,2048+1024,"%02d_%02d_%02d_%02d_%02d_%02d_TLS_w_FEEDBACK_at_Step=%i_of%i_OmL_%.3f_Gamma%.3f_dt_%.4f_CutOff_%.4f.dat",year,month+1,day,hour,minute,second,Nfb,(int)ttotal,Omega,Gamma,dt,cutoff*1000000.);
    file = fopen(FILE_NAME,"w");
    fprintf(file,"#time \tS+S- \t norm \n");
    fprintf(file,"# %02d_%02d_%02d_%02d_%02d_%02d \n",year,month+1,day,hour,minute,second);
    fprintf(file,"%f \t %e \t %e \n",0.*dt,init_22,0.0000);
    //dimension of system bin
    int d = 2; 
    //Create indices and Tensors, s is system, k photon time bin
    //nur wichtig zum Katalogisieren
    auto timebin = IndexType("timebin");
    auto systembin = IndexType("system");
    // make vector of indices
    // gesamter Konfigurationsraum, also alle Freiheitsgrade vorhanden
    // hier Dimension N Zeitschritte (ttotal) + 1 Systemzustand
    auto bin = std::vector<Index>((int)ttotal+2);
    // zeroth index is always system
    // jeder Unterhilbertraum hat jeweils wieder eine bestimmte
    // Anzahl an Freiheitsgraden. Das Systembin hat Dimension 2
    bin.at(1) = Index("s",2,systembin);
    // jeder Zeitbin (also Photonen pro Zeitschritt) hat eine
    // maximale Anzahl an Photonen, oberes Limit muss so eingestellt
    // werden, dass nichts abgeschnitten wird Nbin
    // one timebin for each timestep
    for(int j = 2; j <= ttotal+1; ++j){ bin.at(j) = Index(nameint("k",j),Nbin,timebin); }
    // jetzt habe ich also N+1 Unterhilbertraeume und davon 
    // 1 mit Dimension 2 (System) und N mit Dimensino Nbin
    auto binlink = std::vector<Index>((int)ttotal+1);
    // one timebin for each timestep
    // anscheinend hat System kein Link ... 
    for(int j = 1; j <= ttotal; ++j) { binlink.at(j) = Index(nameint("l",j),1,Link);   }
    // mit gegebenem Konfigurationsraum, MPS kann initialisiert werden
    // initialize MPS
    // Tensor to initialize MPS
    // ich definiere psi als MPS mit ttotal container, die ich dann
    // mit Tensoren fuelle
    MPS psi=MPS(ttotal+1);
    // ich erstelle MPSTensor, um diesen dann in psi einzusetzen
    // eigentlich eine temp-datei
    ITensor MPSTensor = ITensor();
    // create MPS including the systembin
    // Vektor psi hat ttotal freiheitsgrade, Zahl der Unterhilbertraeume
    for(int j = 1; j<=ttotal+1; j++)
    {
       if(j==1)
       {    MPSTensor = ITensor(bin[j],binlink[j]);
            MPSTensor.set(bin[j](2),binlink[j](1),init_22); // Level besetzt im TLS
            MPSTensor.set(bin[j](1),binlink[j](1),init_11); // Level besetzt im TLS
        }
        //right hand vector
        else
        { // der letzte timebin ist nur nach rueckwaerts verlinkt
          if(j==ttotal+1)
          {
             MPSTensor = ITensor(bin[j],binlink[j-1]);
             //initialize timbin in vacuum ... vakuum (1)
             MPSTensor.set(bin[j](1),binlink[j-1](1),1.);
          }
          else
          {
            // all tensors in between get two links
            // anscheinend kann ich in Itensor beliebige viele Links einsetzen
            // er katalogisiert dies automatisch
            MPSTensor = ITensor(bin[j],binlink[j-1],binlink[j]);
            // initialize timbin in vacuum
            MPSTensor.set(bin[j](1),binlink[j-1](1),binlink[j](1),1.);
          }
        }
   // fill MPS with tensors we need to start at 1 because MPS(0) is not defined
   psi.setA((j),MPSTensor);
   }
   
psi.orthogonalize();   
//psi.position(1); // setze Orthocenter auf 1
if ( isOrtho(psi) ) printf("ja\n");
else printf("Nein!\n");
int oc=orthoCenter(psi); printf("Orthocenter at %i \n",oc);
double norm_psi = norm(psi); printf("Norm(psi)=%.5f\n",norm_psi);
   
    auto U_sys =ITensor(bin[1],prime(bin[1]));
    //pumped TLS 
    U_sys.set(bin[1](1),prime(bin[1](2)),-Cplx_i*dt*Omega);
    U_sys.set(bin[1](2),prime(bin[1](1)),-Cplx_i*dt*Omega);


    auto U_dis = ITensor(bin[1],prime(bin[1]),bin[2],prime(bin[2]));
    //H_int with photon creation/destruction in bin k
    for(int phot=1;phot<Nbin;phot++)
    {            
        U_dis.set(bin[1](2),prime(bin[1](1)),bin[2](phot),prime(bin[2](phot+1)),( 1.)*global_fak*sqrt(dt)*Gamma*sqrt(phot));
        U_dis.set(bin[1](1),prime(bin[1](2)),bin[2](phot+1),prime(bin[2](phot)),(-1.)*global_fak*sqrt(dt)*Gamma*sqrt(phot));
    }//phot
    
    auto temp= ITensor(bin[1],prime(bin[1]));
    temp.set(bin[1](1),prime(bin[1](1)),1.);
    temp.set(bin[1](2),prime(bin[1](2)),1.);

    auto U_1  = U_sys*delta(bin[2],prime(bin[2])) + U_dis;
    auto U_2  = (1./2.)  * mapprime(U_1*prime(U_1),2,1);
    auto U_3  = (1./3.)  * mapprime(U_1*prime(U_2),2,1);
    auto U_4  = (1./4.)  * mapprime(U_1*prime(U_3),2,1);
    auto U_5  = (1./5.)  * mapprime(U_1*prime(U_4),2,1);
    auto U_6  = (1./6.)  * mapprime(U_1*prime(U_5),2,1);
    auto U_7  = (1./7.)  * mapprime(U_1*prime(U_6),2,1);
    auto U_8  = (1./8.)  * mapprime(U_1*prime(U_7),2,1);
    auto U_9  = (1./9.)  * mapprime(U_1*prime(U_8),2,1);
    auto U_10 = (1./10.) * mapprime(U_1*prime(U_9),2,1);
    
    //auto U_evo = delta(bin[2],prime(bin[2]))*temp+U_sys*delta(bin[2],prime(bin[2]))+U_dis+ 0.5*mapprime(U_dis*prime(U_dis),2,1);
    auto U_evo = delta(bin[2],prime(bin[2]))*temp + U_1+U_2+U_3+U_4+U_5+U_6+U_7+U_8+U_9+U_10;
      
    ITensor U,S,V,W,SWAP;
    int j=1;
    
    for(j = 1; j<Nfb;j++){//timeloop
        // after first time step change consistently evolutionoperator, switch timebin to timebin+1
        if(j>1){ U_evo=U_evo*delta(bin[j],bin[j+1])*delta(prime(bin[j]),prime(bin[j+1]));  }
        // create combined system state between future timebin and current time bin
        temp=noprime(U_evo * psi.A(j) * psi.A(j+1));
        // set that U is the futurebin, i.e. we change order system - timebin to timebin-system
        if (j==1) U=ITensor(bin[j+1]); 
        // for j>1 we also have to address the pastbin link to the new timebin, otherwise the 
        // system bin carries all its links with it - which is not canonical mps form where on
        // physical index is accompanied with maximal two non-physical links
        // sum_i|(b,i)>U(|psi(2,i)>|(3,b)>)=sum_i|(b,i)>|temp(s,i)(3,b)>=sum_i|(b,i)>|temp(s)(3,b,i)>
        else U=ITensor(bin[j+1],commonIndex(psi.A(j),psi.A(j-1),Link));
        //svd where U and V share the current link 
        // sum_(i,j)|(b,i)>|(b,i,j)>|psi(j,3)>
        svd(temp,U,S,V, {"LeftIndexName", binlink[j].name(),"Cutoff",cutoff,"Maxm",maxm});
        //set orthocenter to V
        V=V*S;
        //now set U as currentbin, and systembin as the next initial state of psi
        // U(|psi(1)|(3,b)>)=Sum_i|(b,i)>|psi(i,3)>
        psi.setA(j,U);
        psi.setA(j+1,V);
        norm_psi=1.-TLS_norm(psi.A(j+1),systembin);
         
       fprintf(file,"%f\t%e\t%e\n",1.*j*dt,TLS_occupation(psi.A(j+1),systembin),norm_psi);
       if ( (j % ((int)Nfb/10)) == 0) printf("time %.3f of %.3f -- norm=%.10f \n",j*dt,t_end*1.,norm_psi);
    }//end timeloop
    
    
// ---- erzeuge neuen zeitevolutionsoperator fuer den Feedbackteil    
    // allocate new time evolution operator 
    // allocate system evo matrix    
    auto U_sysfb =ITensor(bin[1],prime(bin[1]));
    // we distinguish between purely dissipative and feedback part
    // I apply to the next time bin[Nfb+1]=bin[j+1]
    auto U_disfb = ITensor(bin[1],prime(bin[1]),bin[Nfb+1],prime(bin[Nfb+1]));
    // and the feedback bin is the starting bin[2] 
    auto U_fb    = ITensor(bin[1],prime(bin[1]),bin[2],prime(bin[2]));
        
    //pumped TLS 
    U_sysfb.set(bin[1](1),prime(bin[1](2)),-Cplx_i*dt*Omega);
    U_sysfb.set(bin[1](2),prime(bin[1](1)),-Cplx_i*dt*Omega);
    
    //H_int with photon creation/destruction in bin k
    for(int phot=1;phot<Nbin;phot++)
    {            
        U_disfb.set(bin[1](2),prime(bin[1](1)),bin[Nfb+1](phot  ),prime(bin[Nfb+1](phot+1)),global_fak*sqrt(dt)*Gamma*sqrt(phot));
        U_disfb.set(bin[1](1),prime(bin[1](2)),bin[Nfb+1](phot+1),prime(bin[Nfb+1](phot  )),(-1.)*global_fak*sqrt(dt)*Gamma*sqrt(phot));        
//        U_fb.set(bin[1](2),prime(bin[1](1)),bin[2](phot  ),prime(bin[2](phot+1)),(-1.0*global_fak)*sqrt(dt)*Gamma*exp(-Cplx_i*phi)*sqrt(phot));
//        U_fb.set(bin[1](1),prime(bin[1](2)),bin[2](phot+1),prime(bin[2](phot  )),( 1.0*global_fak)*sqrt(dt)*Gamma*exp(Cplx_i*phi)*sqrt(phot));
        U_fb.set(bin[1](2),prime(bin[1](1)),bin[2](phot  ),prime(bin[2](phot+1)),(-1.0*global_fak)*sqrt(dt)*Gamma*sqrt(phot));
        U_fb.set(bin[1](1),prime(bin[1](2)),bin[2](phot+1),prime(bin[2](phot  )),( 1.0*global_fak)*sqrt(dt)*Gamma*sqrt(phot));
    }//phot

    auto U_int = U_disfb*delta(bin[2    ],prime(bin[2    ])) 
                +U_fb   *delta(bin[Nfb+1],prime(bin[Nfb+1]))    ;
        
    auto tempfb = ITensor(bin[1],prime(bin[1]));
    tempfb.set(bin[1](1),prime(bin[1](1)),1.);
    tempfb.set(bin[1](2),prime(bin[1](2)),1.);

    auto U_fb_1  = U_sysfb*delta(bin[Nfb+1],prime(bin[Nfb+1]))*delta(bin[2],prime(bin[2])) +U_int; 
    auto U_fb_2  = (1./2.)  * mapprime(U_fb_1*prime(U_fb_1),2,1);
    auto U_fb_3  = (1./3.)  * mapprime(U_fb_1*prime(U_fb_2),2,1);
    auto U_fb_4  = (1./4.)  * mapprime(U_fb_1*prime(U_fb_3),2,1);
    auto U_fb_5  = (1./5.)  * mapprime(U_fb_1*prime(U_fb_4),2,1);
    auto U_fb_6  = (1./6.)  * mapprime(U_fb_1*prime(U_fb_5),2,1);
    auto U_fb_7  = (1./7.)  * mapprime(U_fb_1*prime(U_fb_6),2,1);
    auto U_fb_8  = (1./8.)  * mapprime(U_fb_1*prime(U_fb_7),2,1);
    auto U_fb_9  = (1./9.)  * mapprime(U_fb_1*prime(U_fb_8),2,1);
    auto U_fb_10 = (1./10.) * mapprime(U_fb_1*prime(U_fb_9),2,1);
    
    auto U_evofb =   delta(bin[Nfb+1],prime(bin[Nfb+1]))*tempfb*delta(bin[2],prime(bin[2])) + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7  + U_fb_8 + U_fb_9 + U_fb_10 ;
                   
// ----- es gibt kein wohldefiniertes Orthocenter solange ich Berechnungen durchfuehre                   
//if ( isOrtho(psi) ) printf("ja\n");
//else printf("Nein!\n");

//int tFB=3;int tSB=3;int tCB=3;
// das bild im MPS sieht aus [pastbin][feedbackbin][systembin][currentbin]
Index iFB,iSB,iCB,iPB; // index fuer feedback, system, und current bin 
Index iFBPB; // index fuer feedback bin und past bin
// zum zeit berechnen
int ehour = loctime -> tm_hour; int eminute = loctime -> tm_min; int esecond = loctime -> tm_sec;
int eyear = 1900+loctime -> tm_year; int emonth = loctime -> tm_mon; int eday  = loctime -> tm_mday;
double time_before,time_now; time_before =esecond + 60*eminute+60*60*ehour;

for(int i=Nfb;i<ttotal;i++)             
{            
    // zuerst den feedbackbin nach i-1 schaufeln         
    for(int k=1;k<Nfb-1;k++) // immer eine Delayschlaufe zurueck
    {
            // k+(i-Nfb) --> 
            SWAP = psi.A(i-Nfb+k)*psi.A(i-Nfb+k+1);
            if (i-Nfb+k==1) U=ITensor(bin[i-Nfb+k+2]); // keine Linksverlinkung fuer psi.A(1)
            else U=ITensor(bin[i-Nfb+k+2],commonIndex(psi.A(i-Nfb+k),psi.A(i-Nfb+k-1),Link)); 
            svd(SWAP,U,S,V,{"LeftIndexName", binlink[i-Nfb+k].name(),"Cutoff",cutoff,"Maxm",maxm});
            psi.setA(i-Nfb+k,U);
            psi.setA(i-Nfb+k+1,S*V);
    }
     // ------ ueberpruefe ob feedbackbin, systembin und currentbin einen Rang von mindestens 2, hoechstens 3 haben
     // tFB=rank(psi.A(i-1)); tSB=rank(psi.A(i)); tCB=rank(psi.A(i+1)); 
     // ------ der Uebersicht wegen identifizierbare die zugehoerigen Links 
     iFB = findtype(psi.A(i-1),timebin); iFBPB = commonIndex(psi.A(i-1),psi.A(i-2));
     iSB = findtype(psi.A(i),systembin);
     iCB = findtype(psi.A(i+1),timebin);
     // ------ kontrahiere FB,SB, und CB und fuehre eine zweifache Schmidtzerlegung durch
     temp=noprime(U_evofb * psi.A(i-1) * psi.A(i) * psi.A(i+1));
     U=ITensor(iCB,iFB,iFBPB);
     svd(temp,U,S,V, {"RightIndexName", binlink[i].name(),"Cutoff",cutoff,"Maxm",maxm});
     // Systembin nach zum naechsten Zeitschritt i+1
     psi.setA(i+1,V);
     // ---- Berechne den Erwartungswert und schiebe dafuer das Orthocenter zum Systembin
     norm_psi=1.-TLS_norm(V*S,systembin);
     fprintf(file,"%f\t%e\t%e\n",1.*i*dt,TLS_occupation(V*S,systembin),norm_psi);
     // zerlege nun currentbin und feedbackbin und stelle feedbackbin nach links und gebe ihm das Orthocenter
     W=U*S; 
     U=ITensor(iFB,iFBPB);
     svd(W,U,S,V, {"RightIndexName", binlink[i-1].name(),"Cutoff",cutoff,"Maxm",maxm});
     // Feedbackbin nach i-1
     psi.setA(i-1,U*S); 
     // Currentbin nach i
     psi.setA(i,V);
     // ----- zurueck swappen, da fuer viele tau-intervalle das hinzuswappen zu lange dauern wuerde ---
     for(int k=1;k<Nfb-1;k++)
     {      
            SWAP = psi.A(i-1-k)*psi.A(i-k);
            if ( (i-k-1)==1 ) U=ITensor(bin[i-Nfb+2]); // kein zu vererbender Link
            else U=ITensor(bin[i-Nfb+2],commonIndex(psi.A(i-k-2),psi.A(i-k-1),Link));
            svd(SWAP,U,S,V,{"LeftIndexName", binlink[i-k-1].name(),"Cutoff",cutoff,"Maxm",maxm});
            // -- lasse Orthocenter beim vorletzten Bin, der im naechsten Schritt wieder zum Systembin geswappt wird
            if (k==Nfb-2) { psi.setA(i-1-k,U);       psi.setA(i-k,V*S); }
            else          { psi.setA(i-1-k,U*S);       psi.setA(i-k,V); }
     }
     // ----- neuen zeitevolutionsoperator erzeugen, und zwar mit neuem feedback bin und neuen currentbin ---- 
     U_evofb=U_evofb*delta(bin[i+1],bin[i+2])*delta(prime(bin[i+1]),prime(bin[i+2]))*delta(bin[i-Nfb+2],bin[i-Nfb+3])*delta(prime(bin[i-Nfb+2]),prime(bin[i-Nfb+3]));  
     // ----- gebe die Normabweichung aus und ueberpruefe ob MPS in kanonischer Form (maximal Tensor dritter Stufe)
     if ( (i % ((int)ttotal/100)) == 0) 
       { 
         curtime = time(NULL); loctime = localtime (&curtime); 
         ehour = loctime -> tm_hour; eminute = loctime -> tm_min; esecond = loctime -> tm_sec;
         time_now = esecond + 60*eminute+60*60*ehour;  
         printf("Feedback time %.3f of %.3f -- Normdeviation=%.10f -- Time elapsed=%.2f sec. \n",dt*i,t_end*1.,norm_psi,time_now-time_before);
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
