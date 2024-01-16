#include "itensor/all.h"

using namespace itensor;

double population(ITensor A, int Nbin)
{
    Index s=findIndex(A,"sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(2),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).real();
}

double coherence(ITensor A, int Nbin)
{
    Index s=findIndex(A,"sys") ;       
    ITensor SpSm = ITensor(s,prime(s));
    SpSm.set(s(1),prime(s(2)),1.);
    return eltC(dag(prime(A,s))*SpSm*A).imag();
}

int main(int argc, char* argv[])
    { 
        
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
    Real init_EE = input.getReal("init_rhoEE_sys");
    Real init_GG = input.getReal("init_rhoGG_sys");
    Real benchmark_factor = input.getReal("benchmark_factor");
    //printf("RhoEE = %.10f \n",init_EE);
    //printf("RhoGG = %.10f \n",init_GG);
    
    FILE * file;
    file = fopen ("g1_output.dat","w");
    fprintf(file,"#time \tS+S- \t norm \n");
    //dimension of system bin
    int d = 2; 
    int gtotal =2*((int)ttotal);

    // ###################################################################################
    // setup physical degrees of freedom bins from 1 to gtotal+1
    auto bin = std::vector<Index>(gtotal+2);
    bin.at(1) = Index(Nbin,"sys");
    for(int j = 2; j <= gtotal+1; ++j){ bin.at(j) = Index(Nbin,"bath"); }
    // setup links in between the physical degrees of freedom
    auto binlink = std::vector<Index>(gtotal+1);
    for(int j = 1; j <= gtotal; ++j) { binlink.at(j) = Index(1,"link");   }
    // ###################################################################################
    // now initialize an MPS 
    MPS psi   = MPS(gtotal+1);
    ITensor MPSTensor = ITensor();
    
    // first bin loaded with system and its initial conditions
    MPSTensor = ITensor(bin[1],binlink[1]);
    MPSTensor.set(bin[1](2),binlink[1](1),init_EE);
    MPSTensor.set(bin[1](1),binlink[1](1),init_GG);
    psi.setA(1,MPSTensor);
    // bins in between 1 and gototal+1 have links left and right and are in the ground state      
    for(int j = 2; j<=gtotal; j++)
    {
      MPSTensor = ITensor(bin[j],binlink[j-1],binlink[j]);
      MPSTensor.set(bin[j](1),binlink[j-1](1),binlink[j](1),1.);
      psi.setA((j),MPSTensor);
    }
    // last bin has links only to the left
    MPSTensor = ITensor(bin[gtotal+1],binlink[gtotal]);
    MPSTensor.set(bin[gtotal+1](1),binlink[gtotal](1),1.);
    psi.setA((gtotal+1),MPSTensor);
    // now I have an MPS, filled with bin[1] bin[2] ... bin[gtotal+1]
    // bin[1] is the system, and bin[1] is changed throughout the algorithm
    // set orthoCenter
    //psi.position(1);
    psi.orthogonalize();
    //printf("%i \n",orthoCenter(psi));

   // ###################################################################################
   // allocate system evo matrix    
    auto H_sys  = ITensor(bin[1],prime(bin[1]));
    //allocate H_int evo matrix - interaction of system with first timebin
    auto H_dis = ITensor(bin[1],prime(bin[1]),bin[2],prime(bin[2]));
        
    //pumped TLS 
    H_sys.set(bin[1](1),prime(bin[1](2)),-Cplx_i*dt*Omega);
    H_sys.set(bin[1](2),prime(bin[1](1)),-Cplx_i*dt*Omega);
    
    //H_int with photon creation/destruction in bin k
    H_dis.set(bin[1](2),prime(bin[1](1)),bin[2](1),prime(bin[2](2)),-Cplx_i*sqrt(dt)*Gamma);
    H_dis.set(bin[1](1),prime(bin[1](2)),bin[2](2),prime(bin[2](1)),-Cplx_i*sqrt(dt)*Gamma);
    
    //dummy Tensor to create "identity matrix" ITensor bug in multiplying two delta Tensors
    auto temp= ITensor(bin[1],prime(bin[1]));
    temp.set(bin[1](1),prime(bin[1](1)),1.);
    temp.set(bin[1](2),prime(bin[1](2)),1.);
    //PrintData(temp);
    
    auto U_evo = delta(bin[2],prime(bin[2]))*temp+H_sys*delta(bin[2],prime(bin[2]))+H_dis+ 0.5*mapPrime(H_dis*prime(H_dis),2,1);
    //PrintData(U_evo);
    // ###################################################################################

    // ###################################################################################
    ITensor U,S,Vp;
    int j;

    double rhoEE=init_EE;
    double rhoGE=0.;
 
    // variables for the exact solution obtained via Laplace transform
    double tGamma = benchmark_factor*Gamma*Gamma;///(2.*3.141);
    double tOmega = 2.*Omega;
    double bOm =sqrt(1.*tOmega*tOmega-tGamma*tGamma/16.);
    double bGam=3.*tGamma/4.;    
    double exact;
    double eG;
    double sOm;
    double cOm;
    double S0 =init_EE+init_GG;
    double W0 =init_EE-init_GG;
    printf("dt=%.10f \n",dt);
    fprintf(file,"%.10f  \t %.10f \t %.10f \t %.10f \t %.10f \n",0.,rhoEE,init_EE,0.,1.0000000);
    
    for(j = 1; j<gtotal-10;j++)
    {
        if(j>1){ U_evo=U_evo*delta(bin[j],bin[j+1])*delta(prime(bin[j]),prime(bin[j+1]));  }

        // fuer psi ####################################################
        temp=noPrime(U_evo * psi.A(j) * psi.A(j+1));
        if (j==1) U=ITensor(bin[j+1]); 
        else U=ITensor(bin[j+1],commonIndex(psi.A(j),psi.A(j-1)));
        svd(temp,U,S,Vp);
        Vp=Vp*S;       
        psi.setA(j,U);
        psi.setA(j+1,Vp);
        psi.position(j+1);
        
        // exact solution of the Mollow problem without initial polarization
        sOm  = sin(dt*j*bOm);
        cOm  = cos(dt*j*bOm);
        eG   = exp(-dt*j*bGam);
        
        if (Omega>0.0000001)
        {    
        exact  = eG*sOm*( -S0     *tGamma/bOm
                          -W0*0.25*tGamma/bOm
                          +S0     *(bGam/bOm)*tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma) );
        exact += eG*cOm*(  S0                *tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma) 
                         + W0 ) ;
        exact +=          -S0                *tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma);
        }
        else exact = -S0+2.*init_EE*exp(-tGamma*j*dt); 

        if ( gtotal > 1000 ) // if number of time steps too large, just save every 1000th data points
        {
        if    ((j % (gtotal/1000)) == 0) 
        { fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",j*dt,population(psi(j+1),Nbin),0.5*exact+0.5*S0,coherence(psi(j+1),Nbin),norm(psi)); }
        } 
        else 
        { 
        fprintf(file,"%.10f \t %.10f \t %.10f \t %.10f \t %.10f \n",j*dt,population(psi(j+1),Nbin),0.5*exact+0.5*S0,coherence(psi(j+1),Nbin),norm(psi)); 
        }
        
        if  ( ((j % (gtotal/20)) == 0) || (j>gtotal+20) ) 
        {  printf("%f \t %f---%f \t %f -- j=%i \n",j*dt,population(psi(j+1),Nbin),0.5*exact+0.5*S0,coherence(psi(j+1),Nbin),j);                        }
                   
    }
    // ###################################################################################    
                
/*    printf(" ------------ and back ------------- \n");
        int i;
        for(i=0;i<gtotal-1;i++)
        {    
        //i=1;
        if (i>0) U_evo=U_evo*delta(bin[j-i],bin[j-(i-1)])*delta(prime(bin[j-i]),prime(bin[j-(i-1)]));
        temp=noPrime(dag(U_evo) * psi.A(j-(i+1)) * psi.A(j-i));
        //PrintData(temp);
        if (i>0) U=ITensor(bin[j-i],commonIndex(psi.A(j-i),psi(j-(i-1))));
        else U=ITensor(bin[j-i] );
        svd(temp,U,S,Vp);
        //PrintData(U);
        //PrintData(Vp);
        Vp=Vp*S;       
        psi.setA(j-i,U);
        psi.setA(j-(i+1),Vp);
        psi.position(j-(i+1));
        rhoEE = population(Vp,Nbin);
        if    ((i % (gtotal/1000)) == 0) fprintf(file,"%.10f \t %.10f \t %.10f \n",(j+i)*dt,rhoEE,norm(psi));
        if  ( ( (i % (gtotal/20)) == 0 ) || (i<20) ) { printf("%f \t %f -- j=%i \n",(j+i)*dt,rhoEE,j); }     
        }
 */       
    fclose(file);    
    return 0;
    }
