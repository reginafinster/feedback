#include "itensor/all.h"
#include "header.h"
 
double liouvillian::trace_check_new(MPS& rho, const std::vector<Index>& s)
{
    // initial value for trace is zero
    double trace = 0.0;
     
    // get size of MPS
    int N = rho.N();
     
    // set up MPO
    MPO MPOTrace = MPO(N);
    ITensor MPOTensor;
     
    for (int i=0; i<N; i++)
    {
        MPOTensor = ITensor(s[i],prime(s[i]));
        MPOTensor.set(s[i](1),prime(s[i](1)),1.0);
        MPOTrace.setA((i+1),MPOTensor);
    }
    
    // call resp. position and apply operator
    Args my_args = Args("Cutoff", 1e-30, "Maxm", 5000);
    MPS *res = new MPS;
    *res = MPS(N);
    MPS rho_temp = liouvillian::apply_nonherm_MPO(MPOTrace, rho, *res, my_args);
    delete res;
    
    // contract rho_MPS
    ITensor temp = rho_temp.A(1);
    for (int k = 1; k<N; k++)
    {
        temp *= rho.A(k);
        temp *= rho_temp.A(k+1);
    }
    temp *= rho.A(N);
     
    // read trace. No scaling needed as no operator has been applied.
    trace = sqrt(temp.real());
     
    return trace;
}
