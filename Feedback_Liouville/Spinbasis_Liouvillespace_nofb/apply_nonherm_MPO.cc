#include "itensor/all.h"
#include "header.h"

using namespace itensor;

MPS liouvillian::apply_nonherm_MPO(MPO& L, MPS& psi, MPS& res, Args args)
// Beachte bei der Eingabe: MPS und MPO müssen mit den selben Site-Indizes versehen sein. Bei beiden muss der 
// Zeilenindex ungeprimed und der Spaltenindex geprimed sein. 
// Ist dem nicht so, wird faelschlicherweise mit dem transponierten MPO multipliziert!
{
    // set up empty MPS for result
    int N = psi.N();
    L.position(1);
    psi.position(1);

    // prime site indices of MPS in order to achieve correct multiplication with non-hermitian MPO
    psi.mapprime(0,1,Site);
    
    // declare indices and tensors
    Index si_left;
    Index si_right;
    Index li_psi;
    Index li_L;
    Index li_res;
    ITensor *temp = new ITensor;
    ITensor *U = new ITensor;
    ITensor *S = new ITensor;
    ITensor *V = new ITensor;
    
     // loop through sites. First and last site are special cases
    for (int i = 1; i<N; i++)
    {
        Index si_left = mapprime(commonIndex(psi.A(i), L.A(i), Site),1,0);
        Index si_right = mapprime(commonIndex(psi.A(i+1), L.A(i+1), Site),1,0);

        // contract tensors depending on their position in the chain
        if (i == 1)
        {
            // define additional link indices
            li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
            li_L = commonIndex(L.A(i+1), L.A(i+2), Link);
            // contract Tensors and decompose them again
            *temp = psi.A(i) * L.A(i) * psi.A(i+1) * L.A(i+1);
            *U = ITensor(si_left);
            *V = ITensor(si_right, li_psi, li_L);
        }
        else if (i == (N-1))
        {
            // define additional link indices
            li_res = commonIndex(res.A(i-1), res.A(i), Link);
            // contract Tensors and decompose them again
            *temp = res.A(i) * psi.A(i+1) * L.A(i+1);
            *U = ITensor(si_left, li_res);
            *V = ITensor(si_right);
        }
        else
        {
            // define additional link indices
            li_res = commonIndex(res.A(i-1), res.A(i), Link);
            li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
            li_L = commonIndex(L.A(i+1), L.A(i+2), Link);
            // contract Tensors and decompose them again
            *temp = res.A(i) * psi.A(i+1) * L.A(i+1);
            *U = ITensor(si_left, li_res);
            *V = ITensor(si_right, li_psi, li_L);
        }
        svd(*temp, *U, *S, *V, args);
        *V = (*V) * (*S);
        //factor(temp, U, V, args);
        
        // store Tensors into resulting MPS:
        res.setA(i,*U);
        res.setA(i+1,*V);
    }//EOF site loop
    
    // restore original state of psi
    psi.mapprime(1,0,Site);
    delete U;
    delete V;
    delete S;
    delete temp;
    
    return res;
}

    
    
//     // loop through sites. First and last site are special cases
//     for (int i = 1; i<N; i++)
//     {
//         if (i==1)
//         {
//             // get indices
//             Index si_left = mapprime(commonIndex(psi.A(i), L.A(i), Site),1,0);
//             Index si_right = mapprime(commonIndex(psi.A(i+1), L.A(i+1), Site),1,0);
//             Index li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
//             Index li_L = commonIndex(L.A(i+1), L.A(i+2), Link);
//             
//             // contract Tensors for the first two sites
//             // and decompose them again. This way, afterwards only one link index 
//             // between two sites exists and result is written into returned MPS
//             temp = psi.A(i) * L.A(i) * psi.A(i+1) * L.A(i+1);
//             Print(temp);
//             U = ITensor(si_left);
//             V = ITensor(si_right, li_psi, li_L);
//             factor(temp, U, V);
//             Print(U);
//             Print(V);
//             
//             // store Tensors into resulting MPS:
//             res.setA(i,U);
//             res.setA(i+1,V);
//             Print(res);
//         }
//         else if (i==(N-1))
//         {
//             // get indices
//             Index si_left = mapprime(commonIndex(psi.A(i), L.A(i), Site),1,0);
//             Index si_right = mapprime(commonIndex(psi.A(i+1), L.A(i+1), Site),1,0);
//             li_res = commonIndex(res.A(i-1), res.A(i), Link);
//             
//             // contract Tensors and decompose them again
//             temp = res.A(i) * psi.A(i+1) * L.A(i+1);
//             Print(temp);
//             U = ITensor(si_left, li_res);
//             V = ITensor(si_right);
//             factor(temp, U, V);
//             Print(U);
//             Print(V);
//             
//             // store Tensors into resulting MPS:
//             res.setA(i,U);
//             res.setA(i+1,V);
//             Print(res);
//         }
//         else
//         {
//             si_left = mapprime(commonIndex(psi.A(i), L.A(i), Site),1,0);
//             si_right = mapprime(commonIndex(psi.A(i+1), L.A(i+1), Site),1,0);
//             li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
//             li_L = commonIndex(L.A(i+1), L.A(i+2), Link);
//             
//             li_res = commonIndex(res.A(i-1), res.A(i), Link);
//             
//             
//             // contract Tensors and use left tensor from above svd, it already contains action of the liouvillian on this one site
//             temp = res.A(i) * psi.A(i+1) * L.A(i+1);
//             Print(temp);
//             U = ITensor(si_left, li_res);
//             V = ITensor(si_right, li_psi, li_L);
//             Print(U);
//             Print(V);
//             factor(temp, U, V);
//             Print(U);
//             Print(V);
//             
//             // store Tensors into resulting MPS:
//             res.setA(i,U);
//             res.setA(i+1,V);
//             Print(res);
//             ////////////////////////////////////////////////////////////////////
//         }









