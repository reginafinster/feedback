#include "mollow.h"


namespace itensor{
    
    Mollow::Mollow(std::string infile)
    {
        // params for system and interaction
        InputGroup input = InputGroup(infile,"input");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        dt = input.getReal("time_step");
        t_end = input.getReal("time_end");
        Omega = input.getReal("Omega_TLS");
        Gamma = input.getReal("Gamma");
        Nbin = input.getInt("Dim_of_bath");
        init_EE = input.getReal("init_rhoEE_sys");
        init_GG = input.getReal("init_rhoGG_sys");
        benchmark_factor = input.getReal("benchmark_factor");
        
//         // params for feedback
//         phi = input.getReal("feedbackphase_in_units_of_phi",0);
//         Nfb = input.getInt("feedback_time");
        
        // calculate some params
        ttotal = t_end;
        //ttotal =2*((int)ttotal);
        // ttotal += Nfb;
        args = Args("Cutoff=",cutoff,"Maxm=",maxm);
        phi *= 3.141592654;
        Mollow::setup_dir_and_filename();
    }
    
    std::string Mollow::get_date_as_string()
    {
        time_t t = time(0); 
        struct tm *now = localtime(&t);
        std::stringstream datestring;
        datestring << now->tm_mday << '-' << (now->tm_mon + 1) << '-' << (now->tm_year + 1900);
        return datestring.str();
    }
    
    ITensor Mollow::identity(Index s)
    {
        ITensor my_identity = ITensor(s,prime(s));
        for(int i=1; i<=s.m(); i++)
        {
            my_identity.set(s(i),prime(s)(i),1.0);
        }
        return my_identity;
    }

    double Mollow::occupation(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(2), prime(s(2)), 1.0);
        ITensor temp = dag(prime(psitensor, s)) * (SpSm*psitensor);
        double occ = temp.real();
        return occ;
    }
    
    double Mollow::ground_occupation(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(1), prime(s(1)), 1.0);
        ITensor temp = dag(prime(psitensor, s)) * (SpSm*psitensor);
        double occ = temp.real();
        return occ;
    }
    
    double Mollow::coherence(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor CoOp = ITensor(s, prime(s));
        CoOp.set(s(1), prime(s(2)), 1.0);
        ITensor temp = dag(prime(psitensor, s)) * (CoOp*psitensor);
        std::complex<double> coherence = temp.cplx();
        return coherence.imag();
    }

    double Mollow::norm(ITensor psi)
    {
        return (dag(psi)*psi).real();
    }

    void Mollow::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        std::string date = Mollow::get_date_as_string();
        std::ostringstream oss;
        // oss << "out/TLS_Mollow_" << date.c_str() << "_Om" << Omega << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phase" << (phi/3.141592654) << ".dat";
        oss << "out/TLS_Mollow_" << date.c_str() << "_Om" << Omega << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << ".dat";
        file_identifier = oss.str();
        datastream = fopen(file_identifier.c_str(), "w");
        fprintf(datastream,"#time \texc level \t ground level \t coherence \t norm \n");
        fclose(datastream);
    }
    
    void Mollow::setup_indexbin()
    {
        bin = std::vector<Index>(ttotal+1);
        
        // first entry of vector is the system
        bin.at(0) = Index("s",2,Site);
        
        // rest is filled with indices for timesteps
        for (int i=1;i<(ttotal+1);i++)
        {
            bin.at(i) = Index(nameint("k",i),Nbin,Site);
        }
    }
    
    void Mollow::initialize_MPS_without_feedback()
    {
        // We need one tensor for each timestep and first one is system
        psi = MPS(ttotal+1);
        
        // set arbitrary init_cond in system, timebins are initialized in vacuum
        ITensor MPSTensor = ITensor();
        for (int i=0;i<(ttotal+1);i++)
        {
            MPSTensor = ITensor(bin[i]);
            if(i==0)
            {
                MPSTensor.set(bin[i](1),init_GG);
                MPSTensor.set(bin[i](2),init_EE);
                println("Initialize system with %.3f in ground and %.3f in upper level\n", init_GG, init_EE);
            }
            else
            {
                MPSTensor.set(bin[i](1),1.0);
            }
            psi.setA((i+1),MPSTensor);
        }
        psi.orthogonalize();
    }
    
//     void Mollow::initialize_MPS_with_feedback()
//     {
//         // We need one tensor for each timestep and one for the system
//         psi = MPS(ttotal+1);
//         
//         // system may be in ground or exited state, timebins are initialized in vacuum
//         // system is at position A(Nfb+1) initially
//         ITensor MPSTensor = ITensor();
// 
//         for (int i=1;i<Nfb+1;i++)
//         {
//             MPSTensor = ITensor(bin[i]);
//             MPSTensor.set(bin[i](1),1.0);
//             psi.setA((i),MPSTensor);
//         }
//         
//         MPSTensor = ITensor(bin[0]);
//         if(init_cond=="ground")
//         {
//             MPSTensor.set(bin[0](1),1.0);
//             println("Initialize system with %.3f in ground and %.3f in upper level\n", init_GG, init_EE);
//         }
//         else if(init_cond=="excited")
//         {
//             MPSTensor.set(bin[0](2),1.0);
//             println("Initialize system in excited state\n");
//         }
//         psi.setA((Nfb+1),MPSTensor);
//         
//         for (int i=Nfb+2;i<ttotal+2;i++)
//         {
//             MPSTensor = ITensor(bin[i-1]);
//             MPSTensor.set(bin[i-1](1),1.0);
//             psi.setA((i),MPSTensor);
//         }
//         psi.orthogonalize();
//     }
    
    void Mollow::set_up_MPO_without_feedback()
    {
        Hsys = ITensor(bin[0],prime(bin[0]));
        Hsys.set(bin[0](1), prime(bin[0](2)),-Cplx_i*Omega*dt);
        Hsys.set(bin[0](2), prime(bin[0](1)),-Cplx_i*Omega*dt);
        
        Hint = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));
        for (int i=1; i<Nbin; i++)
        {
            Hint.set(bin[0](2),prime(bin[0](1)),bin[1](i),prime(bin[1](i+1)), -Cplx_i*Gamma*pow(dt,0.5)*pow(i,0.5));
            Hint.set(bin[0](1),prime(bin[0](2)),bin[1](i+1),prime(bin[1](i)), -Cplx_i*Gamma*pow(dt,0.5)*pow(i,0.5));
        }
        
        ITensor identity_system = Mollow::identity(bin[0]);
        ITensor identity_int = Mollow::identity(bin[1]);
        
//         ITensor U_1  = Hsys*identity_int + Hint;
//         ITensor U_2  = (1./2.)  * mapprime(U_1*prime(U_1),2,1);
//         ITensor U_3  = (1./3.)  * mapprime(U_1*prime(U_2),2,1);
//         ITensor U_4  = (1./4.)  * mapprime(U_1*prime(U_3),2,1);
//         ITensor U_5  = (1./5.)  * mapprime(U_1*prime(U_4),2,1);
//         ITensor U_6  = (1./6.)  * mapprime(U_1*prime(U_5),2,1);
//         ITensor U_7  = (1./7.)  * mapprime(U_1*prime(U_6),2,1);
//         ITensor U_8  = (1./8.)  * mapprime(U_1*prime(U_7),2,1);
//         ITensor U_9  = (1./9.)  * mapprime(U_1*prime(U_8),2,1);
//         ITensor U_10 = (1./10.) * mapprime(U_1*prime(U_9),2,1);
//         
//         Uevo = identity_int*identity_system + U_1+U_2+U_3+U_4+U_5+U_6+U_7+U_8+U_9+U_10;
        
        Uevo = identity_int*identity_system + Hsys*identity_int + Hint;
        //Uevo = identity_int*identity_system + Hsys*identity_int + Hint + 0.5*mapprime(prime(Hint)*Hint,2,1);
    }
    
//     void Mollow::set_up_MPO_with_feedback()
//     {
//         Hsys = ITensor(bin[0],prime(bin[0]));
//         Hsys.set(bin[0](1), prime(bin[0](2)),-Cplx_i*Omega*dt);
//         Hsys.set(bin[0](2), prime(bin[0](1)),-Cplx_i*Omega*dt);
//         
//         Hint = ITensor(bin[0],prime(bin[0]),bin[Nfb+1],prime(bin[Nfb+1]));
//         Hfb = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));
// 
//         for(int i=1;i<Nbin;i++)
//         {            
//             Hint.set(bin[0](2),prime(bin[0](1)),bin[Nfb+1](i),prime(bin[Nfb+1](i+1)), pow(dt,0.5)*Gamma*pow(i,0.5));
//             Hint.set(bin[0](1),prime(bin[0](2)),bin[Nfb+1](i+1),prime(bin[Nfb+1](i)), (-1.)*pow(dt,0.5)*Gamma*pow(i,0.5));        
//             Hfb.set(bin[0](2),prime(bin[0](1)),bin[1](i),prime(bin[1](i+1)), (-1.0)*pow(dt,0.5)*Gamma*exp(-Cplx_i*phi)*pow(i,0.5));
//             Hfb.set(bin[0](1),prime(bin[0](2)),bin[1](i+1),prime(bin[1](i)), pow(dt,0.5)*Gamma*exp(Cplx_i*phi)*pow(i,0.5));
//         }
//         
//         ITensor identity_system = Mollow::identity(bin[0]);
//         ITensor identity_int = Mollow::identity(bin[Nfb+1]);
//         ITensor identity_fb = Mollow::identity(bin[1]);
//         
//         Hsys *= (identity_int*identity_fb);
//         Hint *= identity_fb;
//         Hfb *= identity_int;
//         
//         Hint += Hfb;
//         
//         ITensor U_fb_1  = Hsys + Hint; 
//         ITensor U_fb_2  = (1./2.)  * mapprime(U_fb_1*prime(U_fb_1),2,1);
//         ITensor U_fb_3  = (1./3.)  * mapprime(U_fb_1*prime(U_fb_2),2,1);
//         ITensor U_fb_4  = (1./4.)  * mapprime(U_fb_1*prime(U_fb_3),2,1);
//         ITensor U_fb_5  = (1./5.)  * mapprime(U_fb_1*prime(U_fb_4),2,1);
//         ITensor U_fb_6  = (1./6.)  * mapprime(U_fb_1*prime(U_fb_5),2,1);
//         ITensor U_fb_7  = (1./7.)  * mapprime(U_fb_1*prime(U_fb_6),2,1);
//         ITensor U_fb_8  = (1./8.)  * mapprime(U_fb_1*prime(U_fb_7),2,1);
//         ITensor U_fb_9  = (1./9.)  * mapprime(U_fb_1*prime(U_fb_8),2,1);
//         ITensor U_fb_10 = (1./10.) * mapprime(U_fb_1*prime(U_fb_9),2,1);
//     
//         Uevo = identity_system*identity_int*identity_fb + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7 + U_fb_8 + U_fb_9 + U_fb_10 ;
//         
//         // second order evolution
//         // Hint += Hfb;
//         // Uevo = identity_system*identity_int*identity_fb + Hsys + Hint + 0.5*mapprime(prime(Hint)*Hint,2,1);
//     }
    
    void Mollow::timeevolve_system_without_feedback()
    {
        datastream = fopen(file_identifier.c_str(), "a");
        Mollow::setup_indexbin();
        Mollow::initialize_MPS_without_feedback();
        Mollow::set_up_MPO_without_feedback();
        
        // set up variables for exact solution
        tGamma = benchmark_factor*Gamma*Gamma;
        tOmega = 2.0*Omega;
        bOm = pow(tOmega*tOmega-tGamma*tGamma/16.0, 0.5);
        bGam=3.0*tGamma/4.0;
        S0 =init_EE+init_GG;
        W0 =init_EE-init_GG;
        
        // time evolve system
        for(int j = 0; j<ttotal; j++)
        {
            Mollow::apply_Uevo_without_feedback(j);
            
            if(j<(ttotal-1))
            {
                Uevo = Uevo*delta(bin[j+1],bin[j+2])*delta(prime(bin[j+1]),prime(bin[j+2]));
            }
            
            Mollow::calc_exact_variables(j);
            
            if ((j%10)==0) 
            {
                printf("time %.3f of %.3f -- norm=%.10f \n", j*dt, t_end*1., Mollow::norm(psi.A(j+2)));
                fprintf(datastream, "%f \t %e \t %f \t %e \t %e \n", j*dt, Mollow::occupation(psi.A(j+2), Site), Mollow::ground_occupation(psi.A(j+2), Site), Mollow::coherence(psi.A(j+2),Site), Mollow::norm(psi.A(j+2)));
            }
        }
        
        fclose(datastream);
    }
    
    void Mollow::apply_Uevo_without_feedback(int timestep)
    {
        ITensor U,S,V,temp;    
        temp = noprime(Uevo*psi.A(timestep+1)*psi.A(timestep+2));
        U = ITensor(bin[timestep+1],commonIndex(psi.A(timestep+1),psi.A(timestep),Link));
        V = ITensor(bin[0],commonIndex(psi.A(timestep+2),psi.A(timestep+3)));
        svd(temp,U,S,V);
        psi.setA(timestep+1,U);
        psi.setA(timestep+2,V*S);
    }
    
    void Mollow::calc_exact_variables(int timestep)
    {
        // exact solution of the Mollow problem without initial polarization
        sOm  = sin(dt*timestep*bOm);
        cOm  = cos(dt*timestep*bOm);
        eG   = exp(-dt*timestep*bGam);
        
        if (Omega>0.0000001)
        {    
        exact  = eG*sOm*( -S0     *tGamma/bOm
                        -W0*0.25*tGamma/bOm
                        +S0     *(bGam/bOm)*tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma) );
        exact += eG*cOm*(  S0                *tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma) 
                        + W0 ) ;
        exact +=          -S0                *tGamma*tGamma/(2.*tOmega*tOmega+tGamma*tGamma);
        }
        else exact = -S0+2.*init_EE*exp(-tGamma*timestep*dt);
    }
    
//     void Mollow::timeevolve_system_with_feedback()
//     {
//         datastream = fopen(file_identifier.c_str(), "a");
//         
//         Mollow::setup_indexbin();
//         Mollow::initialize_MPS();
//         printf("time %.3f of %.3f -- norm=%.10f -- occ=%.10f \n", 0.0*dt, t_end, (1.0 - Mollow::norm(psi.A(1), Site)), Mollow::TLS_occupation(psi.A(1), Site));
//         Mollow::set_up_MPO_with_feedback();
//         for(int j=Nfb; j<ttotal; j++)
//         {
//             Mollow::apply_Uevo_with_feedback(j);
//             if(j<(ttotal-1))
//             {
//                 // Uevo needs to act on future bin of MPS for next timestep:
//                 Uevo *= delta(bin[j+1],bin[j+2]);
//                 Uevo *= delta(prime(bin[j+1]),prime(bin[j+2]));
//                 Uevo *= delta(bin[j-Nfb+1],bin[j-Nfb+2]);
//                 Uevo *= delta(prime(bin[j-Nfb+1]),prime(bin[j-Nfb+2]));
//             }
//         }
//         fclose(datastream);
//     }
    
//     void Mollow::apply_Uevo_with_feedback(int timestep)
//     {
//         Mollow::swap_feedbackbin_to_system(timestep);
//         
//         ITensor U,S,V,temp;
//         //shift ortho-center from feedbackbin to systembin to compute observable
//         temp = psi.A(timestep)*psi.A(timestep+1);
//         U = ITensor(bin[0],rightLinkInd(psi,timestep+1));
//         svd(temp,U,S,V,args);
//         psi.setA(timestep,V);
//         psi.setA(timestep+1,U*S);
//         fprintf(datastream, "%f \t %e \t %e \n", (timestep-Nfb)*dt, Mollow::TLS_occupation(psi.A(timestep+1), Site), (1.0 - Mollow::norm(psi.A(timestep+1), Site)));
//         if ((timestep%10)==0) 
//         {
//             printf("time %.3f of %.3f -- norm=%.10f \n",(timestep-Nfb)*dt,t_end*1., (1.0 - Mollow::norm(psi.A(timestep+1), Site)));
//         }
//         
//         Index currenttimebin = findtype(psi.A(timestep+2),Site);
//         
//         temp = noprime(psi.A(timestep)*psi.A(timestep+1)*psi.A(timestep+2)*Uevo);
//         
//         // bring system to A(timestep+2), thus to the left
//         U = ITensor(bin[0],leftLinkInd(psi,timestep+3));
//         svd(temp,U,S,V,args);
//         V *= S;
//         psi.setA(timestep+2,U);
//             
//         // bring current timebin one position to the left (which is middle position)
//         // and feedbacktimebin back to its position at A(timestep) 
//         U = ITensor(currenttimebin,commonIndex(U,V));
//         temp=V;
//         svd(temp,U,S,V,args);
//         psi.setA(timestep+1,U);
//         psi.setA(timestep,V*S);
//         
//         Mollow::swap_feedbackbin_to_origin(timestep);
//     }
//     
//     void Mollow::swap_feedbackbin_to_system(int timestep)
//     {
//         ITensor U,S,V,temp;
//         Index feedbacktimebin = findtype(psi.A(timestep+1-Nfb),Site);
//         for(int i=timestep+1-Nfb; i<timestep; i++)
//         {
//             temp = psi.A(i)*psi.A(i+1);
//             U = ITensor(feedbacktimebin, commonIndex(psi.A(i+1),psi.A(i+2)));
//             svd(temp,U,S,V,args);
//             psi.setA(i,V);
//             psi.setA(i+1,U*S);
//         }
//     }
//     
//     void Mollow::swap_feedbackbin_to_origin(int timestep)
//     {
//         ITensor temp,U,S,V;
//         Index feedbacktimebin = findtype(psi.A(timestep),Site);
//         // loop ends when feedback bin is located one bin right of its original position
//         // because in the last step, the orthocenter needs to be placed
//         // at the future feedback bin for the next time step
//         for(int i=timestep; i>timestep-Nfb+2; i--) 
//         {
//             temp = psi.A(i-1)*psi.A(i);
//             U = ITensor(feedbacktimebin, leftLinkInd(psi,i-1));
//             svd(temp,U,S,V,args);
//             psi.setA(i,V);
//             psi.setA(i-1,U*S);
//         }
//         temp = psi.A(timestep-Nfb+1)*psi.A(timestep-Nfb+2); 
//         U = ITensor(feedbacktimebin, leftLinkInd(psi,timestep-Nfb+1));
//         svd(temp,U,S,V,args);
//         psi.setA(timestep-Nfb+2,V*S);
//         psi.setA(timestep-Nfb+1,U);
//     }

}//eof namespace




