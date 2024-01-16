#include "TLS_Liouville.h"


namespace itensor{
    
    TLS_Liouville::TLS_Liouville(std::string infile)
    {
        // params for time evolution
        InputGroup input = InputGroup(infile,"input");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        dt = input.getReal("time_step");
        t_end = input.getReal("time_end");
        Omega = input.getReal("Omega_TLS");
        init_cond = input.getString("initial_condition");
        args = Args("Cutoff=",cutoff,"Maxm=",maxm);
        
        // params for interaction
        Gamma = input.getReal("Gamma");
        LindbladLossGamma = input.getReal("LindbladLossGamma");
        LindbladPumpGamma = input.getReal("LindbladPumpGamma");
        GammaPureDephasing = input.getReal("GammaPureDephasing");
        
        // params for feedback
        Nfb = input.getInt("feedback_time");
        phi = input.getReal("feedbackphase_in_units_of_phi",0.0);
        
        // params for spectrum
        spectrum_steps = input.getInt("spectrum_steps");
        spectrum_interval = input.getReal("spectrum_interval");  
        
        // params for program control
        save_nth_timestep = input.getInt("save_nth_timestep");
        print_nth_timestep = input.getInt("print_nth_timestep");
        use_feedback = input.getString("use_feedback");
        
        // calculate some params
        ttotal = t_end + Nfb;
        phi *= 3.141592654;
        
        LindbladLossGamma /= 3.141592654;
        LindbladLossGamma = pow(LindbladLossGamma,2.0);
        LindbladLossGamma *= 2.0;
        
        LindbladPumpGamma /= 3.141592654;
        LindbladPumpGamma = pow(LindbladPumpGamma,2.0);
        LindbladPumpGamma *= 2.0;
        
        if (spectrum_steps >= t_end)
        {
            println("spectrum_steps must be smaller than number of timesteps. Exiting");
            throw std::exception();
        }

        TLS_Liouville::setup_dir_and_filename();
    }
    
    std::string TLS_Liouville::get_date_as_string()
    {
        time_t t = time(0); 
        struct tm *now = localtime(&t);
        std::stringstream datestring;
        datestring << now->tm_mday << '-' << (now->tm_mon + 1) << '-' << (now->tm_year + 1900);
        return datestring.str();
    }
    
    ITensor TLS_Liouville::identity(Index s)
    {
        ITensor my_identity = ITensor(s,prime(s));
        for(int i=1; i<=s.m(); i++)
        {
            my_identity.set(s(i),prime(s)(i),1.0);
        }
        return my_identity;
    }

    double TLS_Liouville::excited_TLS_occupation(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(4), prime(s(4)), 1.0);
        ITensor temp = dag(prime(psitensor, Site)) * (SpSm*psitensor);
        double occ = temp.real();
        occ = pow(occ,0.5);
        return occ;
    }
    
    double TLS_Liouville::ground_TLS_occupation(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(1), prime(s(1)), 1.0);
        ITensor temp = dag(prime(psitensor, Site)) * (SpSm*psitensor);
        double occ = temp.real();
        occ = pow(occ,0.5);
        return occ;
    }
    
    double TLS_Liouville::trace(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor TraceOp1 = ITensor(s, prime(s));
        TraceOp1.set(s(1), prime(s(1)), 1.0);
        ITensor temp1 = dag(prime(psitensor, Site)) * (TraceOp1*psitensor);
        double trace1 = temp1.real();
        trace1 = pow(trace1,0.5);
        
        ITensor TraceOp2 = ITensor(s, prime(s));
        TraceOp2.set(s(4), prime(s(4)), 1.0);
        ITensor temp2 = dag(prime(psitensor, Site)) * (TraceOp2*psitensor);
        double trace2 = temp2.real();
        trace2 = pow(trace2,0.5);
        
        double trace = trace1 + trace2;
        return trace;
    }
    
    ITensor contract_psi(MPS psi, int ttotal)
    {
        ITensor temp = psi.A(1);
        for (int j = 2; j<=(ttotal+1); j++)
        {
            temp *= psi.A(j);
        }
        return temp;
    
    }

    double TLS_Liouville::rhosquare(ITensor psi, IndexType type)
    {
        return (dag(psi)*psi).real();
    }
    
    ITensor TLS_Liouville::annihilator(Index s)
    {
        ITensor Sm = ITensor(s, prime(s));
        Sm.set(s(2), prime(s(2)), 1.0);
        return Sm;
    }
    
    std::complex<double> TLS_Liouville::calc_g1_bath(Index timebin, ITensor temp)
    {
        ITensor Sm_temp = TLS_Liouville::annihilator(timebin);
        g1_tensor = dag(noprime(temp*Sm_RefBin)) * noprime(Sm_temp*temp);
        std::complex<double> g1_cplx = g1_tensor.cplx();
        return g1_cplx;
    }
    
    std::complex<double> TLS_Liouville::calc_g2_bath(Index timebin, ITensor temp)
    {
        ITensor Sm_temp = TLS_Liouville::annihilator(timebin);
        g2_tensor = dag(noprime(Sm_temp*noprime(Sm_RefBin*temp))) * noprime(Sm_temp*noprime(Sm_RefBin*temp));
        std::complex<double> g2_cplx = g2_tensor.cplx();
        return g2_cplx;
    }
    
    std::complex<double> TLS_Liouville::calc_pop_bath(ITensor temp)
    {
        pop_bath_tensor = dag(noprime(Sm_RefBin*temp)) * noprime(Sm_RefBin*temp);            
        std::complex<double> pop_bath_cplx = pop_bath_tensor.cplx();
        return pop_bath_cplx;
    }
    
    void TLS_Liouville::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::string date = TLS_Liouville::get_date_as_string();
        std::ostringstream oss1, oss2, oss3;
        oss1 << "out/timeevo_" << date.c_str() << "fb_" << use_feedback << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamL"<< LindbladLossGamma << "_GamP" << LindbladPumpGamma << "_GamD" << GammaPureDephasing << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phi" << (phi/3.141592654) << ".dat";
        file_identifier_timeevo = oss1.str();
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "w");
        fprintf(datastream_timeevo,"#time \tS+S- \tS-S+ \t trace \t purity \n");
        fclose(datastream_timeevo);
        
        oss2 << "out/corr_" << date.c_str() << "fb_" << use_feedback << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamL"<< LindbladLossGamma << "_GamP" << LindbladPumpGamma << "_GamD" << GammaPureDephasing << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phi" << (phi/3.141592654) << ".dat";
        file_identifier_correlation = oss2.str();
        datastream_correlation = fopen(file_identifier_correlation.c_str(), "w");
        fprintf(datastream_correlation, "#time \t Re[g1] \t Im[g1] \t Re[g2] \t Im[g2] \t pop_bath^2 \n");
        fclose(datastream_correlation);
        
        oss3 << "out/spec_" << date.c_str() << "fb_" << use_feedback << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamL"<< LindbladLossGamma << "_GamP" << LindbladPumpGamma << "_GamD" << GammaPureDephasing << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phi" << (phi/3.141592654) << ".dat";
        file_identifier_spectrum = oss3.str();
        datastream_spectrum = fopen(file_identifier_spectrum.c_str(), "w");
        fprintf(datastream_spectrum, "#omega \t Re[S] \t Im[S] \t bath_pop \n");
        fclose(datastream_spectrum);
    }
    
    void TLS_Liouville::setup_indexbin()
    {
        bin = std::vector<Index>(ttotal+1);
        
        // first entry of vector is the system
        bin.at(0) = Index("s", 4, Site);
        
        // rest is filled with indices for timesteps
        for (int i=1;i<(ttotal+1);i++)
        {
            bin.at(i) = Index(nameint("k",i),4,Site);
        }
    }
    
    void TLS_Liouville::initialize_MPS_with_feedback()
    {
        // We need one tensor for each timestep and first one is system
        psi = MPS(ttotal+1);
        
        // system may be in ground or exited state, timebins are initialized in vacuum
        // system is at position A(Nfb+1) initially
        ITensor MPSTensor = ITensor();

        for (int i=1;i<Nfb+1;i++)
        {
            MPSTensor = ITensor(bin[i]);
            MPSTensor.set(bin[i](1),1.0);
            psi.setA((i),MPSTensor);
        }
        
        MPSTensor = ITensor(bin[0]);
        if(init_cond=="ground")
        {
            MPSTensor.set(bin[0](1),1.0);
            MPSTensor.set(bin[0](2),0.5);
            println("Initialize system in ground state\n");
        }
        else if(init_cond=="excited")
        {
            MPSTensor.set(bin[0](4),1.0);
            println("Initialize system in excited state\n");
        }
        psi.setA((Nfb+1),MPSTensor);
        
        for (int i=Nfb+2;i<ttotal+2;i++)
        {
            MPSTensor = ITensor(bin[i-1]);
            MPSTensor.set(bin[i-1](1),1.0);
            psi.setA((i),MPSTensor);
        }
        // orthogonalize in order to establish links
        // and shift orthocenter to system when no feedback is used
        if (use_feedback=="yes")
        {
            psi.orthogonalize();
        }
        else if (use_feedback=="no")
        {
            psi.orthogonalize();
            psi.position(Nfb+1);
        }
    }
    
    void TLS_Liouville::set_up_Uevo_with_feedback()
    {
        Lsys = ITensor(bin[0],prime(bin[0]),bin[Nfb+1],prime(bin[Nfb+1]));
        for (int p=1; p<5; p++)
        {
            Lsys.set(bin[0](1),bin[Nfb+1](p),prime(bin[0](3)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](2),bin[Nfb+1](p),prime(bin[0](4)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](3),bin[Nfb+1](p),prime(bin[0](1)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](4),bin[Nfb+1](p),prime(bin[0](2)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](2),bin[Nfb+1](p),prime(bin[0](1)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](4),bin[Nfb+1](p),prime(bin[0](3)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](1),bin[Nfb+1](p),prime(bin[0](2)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](3),bin[Nfb+1](p),prime(bin[0](4)),prime(bin[Nfb+1](p)),-Cplx_i*Omega*dt);
        }
        
        Lint = ITensor(bin[0],prime(bin[0]),bin[Nfb+1],prime(bin[Nfb+1]));
        for (int p=1; p<=2; p++)
        {
            Lint.set(bin[0](1),bin[Nfb+1](p+2),prime(bin[0](3)),prime(bin[Nfb+1](p)), Gamma*pow(dt,0.5));
            Lint.set(bin[0](2),bin[Nfb+1](p+2),prime(bin[0](4)),prime(bin[Nfb+1](p)), Gamma*pow(dt,0.5));
        }
        for (int p=3; p<=4; p++)
        {
            Lint.set(bin[0](3),bin[Nfb+1](p-2),prime(bin[0](1)),prime(bin[Nfb+1](p)), -Gamma*pow(dt,0.5));
            Lint.set(bin[0](4),bin[Nfb+1](p-2),prime(bin[0](2)),prime(bin[Nfb+1](p)), -Gamma*pow(dt,0.5));
        }
        for (int p=2; p<=4; p+=2)
        {
            Lint.set(bin[0](2),bin[Nfb+1](p-1),prime(bin[0](1)),prime(bin[Nfb+1](p)), -Gamma*pow(dt,0.5));
            Lint.set(bin[0](4),bin[Nfb+1](p-1),prime(bin[0](3)),prime(bin[Nfb+1](p)), -Gamma*pow(dt,0.5));
        }
        for (int p=1; p<=3; p+=2)
        {
            Lint.set(bin[0](1),bin[Nfb+1](p+1),prime(bin[0](2)),prime(bin[Nfb+1](p)), Gamma*pow(dt,0.5));
            Lint.set(bin[0](3),bin[Nfb+1](p+1),prime(bin[0](4)),prime(bin[Nfb+1](p)), Gamma*pow(dt,0.5));
        }

        Lfb = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));
        for (int p=1; p<=2; p++)
        {
            Lfb.set(bin[0](1),bin[1](p+2),prime(bin[0](3)),prime(bin[1](p)), -Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](2),bin[1](p+2),prime(bin[0](4)),prime(bin[1](p)), -Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
        }
        for (int p=3; p<=4; p++)
        {
            Lfb.set(bin[0](3),bin[1](p-2),prime(bin[0](1)),prime(bin[1](p)), Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](4),bin[1](p-2),prime(bin[0](2)),prime(bin[1](p)), Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
        }
        for (int p=2; p<=4; p+=2)
        {
            Lfb.set(bin[0](2),bin[1](p-1),prime(bin[0](1)),prime(bin[1](p)), Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](4),bin[1](p-1),prime(bin[0](3)),prime(bin[1](p)), Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
        }
        for (int p=1; p<=3; p+=2)
        {
            Lfb.set(bin[0](1),bin[1](p+1),prime(bin[0](2)),prime(bin[1](p)), -Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](3),bin[1](p+1),prime(bin[0](4)),prime(bin[1](p)), -Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
        }
        
        LindbladLoss = ITensor(bin[0],prime(bin[0]));
        LindbladLoss.set(bin[0](4),prime(bin[0](1)), LindbladLossGamma);
        LindbladLoss.set(bin[0](2),prime(bin[0](2)), -0.5*LindbladLossGamma);
        LindbladLoss.set(bin[0](3),prime(bin[0](3)), -0.5*LindbladLossGamma);
        LindbladLoss.set(bin[0](4),prime(bin[0](4)), -1.0*LindbladLossGamma);
        
        LindbladPump = ITensor(bin[0],prime(bin[0]));
        LindbladPump.set(bin[0](1),prime(bin[0](4)), LindbladPumpGamma);
        LindbladPump.set(bin[0](1),prime(bin[0](1)), -1.0*LindbladPumpGamma);
        LindbladPump.set(bin[0](2),prime(bin[0](2)), -0.5*LindbladPumpGamma);
        LindbladPump.set(bin[0](3),prime(bin[0](3)), -0.5*LindbladPumpGamma);
        
        Lindblad = ITensor(bin[0],prime(bin[0]));
        Lindblad = LindbladLoss + LindbladPump;
        
        ITensor identity_system = TLS_Liouville::identity(bin[0]);
        ITensor identity_int = TLS_Liouville::identity(bin[Nfb+1]);
        ITensor identity_fb = TLS_Liouville::identity(bin[1]);
        
        Lindblad *= (identity_int*identity_fb);
        Lsys *= identity_fb;
        Lfb *= identity_int;
        Lint *= identity_fb;
        
        ITensor U_fb_1  = Lsys + Lint + Lfb + Lindblad; 
        ITensor U_fb_2  = (1./2.)  * mapprime(U_fb_1*prime(U_fb_1),2,1);
        ITensor U_fb_3  = (1./3.)  * mapprime(U_fb_1*prime(U_fb_2),2,1);
        ITensor U_fb_4  = (1./4.)  * mapprime(U_fb_1*prime(U_fb_3),2,1);
        ITensor U_fb_5  = (1./5.)  * mapprime(U_fb_1*prime(U_fb_4),2,1);
        ITensor U_fb_6  = (1./6.)  * mapprime(U_fb_1*prime(U_fb_5),2,1);
        ITensor U_fb_7  = (1./7.)  * mapprime(U_fb_1*prime(U_fb_6),2,1);
        ITensor U_fb_8  = (1./8.)  * mapprime(U_fb_1*prime(U_fb_7),2,1);
        ITensor U_fb_9  = (1./9.)  * mapprime(U_fb_1*prime(U_fb_8),2,1);
        ITensor U_fb_10 = (1./10.) * mapprime(U_fb_1*prime(U_fb_9),2,1);
    
        Uevo = identity_system*identity_int*identity_fb + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7 + U_fb_8 + U_fb_9 + U_fb_10;

        // second order evolution
        // Levo = identity_int*identity_system + Lsys + Lint + 0.5*mapprime(prime(Lint)*Lint,2,1);
    }
    
    void TLS_Liouville::apply_Uevo_with_feedback(int timestep)
    {
        // wenn kein Feedback: orthocenter auf system
        // wenn Feedback: orthocenter auf feedbackbin
        ITensor U,S,V,temp;
        if (use_feedback=="yes")
        {
            // orthocenter is at feedbackbin from previous timestep
            // bring feedbackbin to the left of the systembin
            // thus to A(timestep)
            TLS_Liouville::swap_feedbackbin_to_system(timestep);
            //shift ortho-center from feedbackbin to systembin 
            temp = psi.A(timestep)*psi.A(timestep+1);
            U = ITensor(bin[0],rightLinkInd(psi,timestep+1));
            svd(temp,U,S,V,args);
            psi.setA(timestep,V);
            psi.setA(timestep+1,U*S);
            // system is now at A(timestep+1)
        }
        
        // compute observables
        
        // contract MPS and read out trace and occupation
        ITensor psi_contracted = TLS_Liouville::contract_psi(psi,ttotal);
        double sys_ground_occ = psi_contracted.real(s[0](1),s[1] );
        
        if ((timestep%save_nth_timestep)==0)
        {
            fprintf(datastream_timeevo, "%f \t %e \t %e \t %e \t %e \n", (timestep-Nfb)*dt, TLS_Liouville::excited_TLS_occupation(psi.A(timestep+1), Site), TLS_Liouville::ground_TLS_occupation(psi.A(timestep+1), Site), TLS_Liouville::trace(psi.A(timestep+1), Site), TLS_Liouville::rhosquare(psi.A(timestep+1), Site));
        }
        if ((timestep%print_nth_timestep)==0) 
        {
            printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", (timestep-Nfb)*dt, t_end*dt, TLS_Liouville::excited_TLS_occupation(psi.A(timestep+1), Site),  TLS_Liouville::trace(psi.A(timestep+1), Site) ,TLS_Liouville::rhosquare(psi.A(timestep+1), Site));
        }
        
        // bins are now located: feedbackbin - system - currentbin

        if (use_feedback=="yes")
        {
            // apply Uevo:
            temp = noprime(psi.A(timestep)*psi.A(timestep+1)*psi.A(timestep+2)*Uevo);
            
            // temp has three site-indices, thus we need to decompose twice
            // system is brought to A(timestep+2), one to the right
            U = ITensor(bin[0],leftLinkInd(psi,timestep+3));
            svd(temp,U,S,V,args);
            V *= S;
            psi.setA(timestep+2,U);
                
            // bring current timebin one position to the left (which is middle position)
            // and feedbacktimebin back to its position at A(timestep) 
            Index currenttimebin = findtype(psi.A(timestep+2),Site);
            U = ITensor(currenttimebin,commonIndex(U,V));
            temp=V;
            svd(temp,U,S,V,args);
            psi.setA(timestep+1,U);
            psi.setA(timestep,V*S);
            
            TLS_Liouville::swap_feedbackbin_to_origin(timestep);
        }
        else if (use_feedback=="no")
        {
            // create dummy vaccum timebin 
            Index dummy_vacuum_site_index = findtype(psi.A(timestep+1-Nfb),Site);
            Index dummylink1 = leftLinkInd(psi,timestep+1-Nfb);
            Index dummylink2 = rightLinkInd(psi,timestep+1-Nfb);
            ITensor dummy_vacuum = ITensor(dummy_vacuum_site_index,dummylink1,dummylink2);
            dummy_vacuum.set(dummy_vacuum_site_index(1), dummylink1(1),dummylink2(1),1.0);
            
            // apply Uevo:
            temp = noprime(dummy_vacuum*psi.A(timestep+1)*psi.A(timestep+2)*Uevo);
            
            // temp has three site-indices, thus we need to decompose twice
            // system is brought to A(timestep+2), one to the right
            U = ITensor(bin[0],leftLinkInd(psi,timestep+3));
            svd(temp,U,S,V,args);
            V *= S;
            psi.setA(timestep+2,U);

            // bring current timebin one position to the left of system (which is middle position at A(timestep+1)
            // and discard the dummy vaccum bin
            Index currenttimebin = findtype(psi.A(timestep+2),Site);
            U = ITensor(currenttimebin,commonIndex(U,V));
            temp=V;
            svd(temp,U,S,V,args);
            psi.setA(timestep+1,U*S);
            
            // shift orthocenter to system for next timestep
            // thus one to the right to A(timestep+2)
            temp = psi.A(timestep+1)*psi.A(timestep+2);
            U = ITensor(bin[0],rightLinkInd(psi,timestep+2));
            svd(temp,U,S,V,args);
            psi.setA(timestep+2,U*S);
            psi.setA(timestep+1,V);
        }
        else
        {
             println("Could not identify use of feedback. Exiting");
             throw std::exception();
        }
        // positions are : feedbackbin - currentbin - systembin with orthocenter
    }
    
    void TLS_Liouville::swap_feedbackbin_to_system(int timestep)
    {
        ITensor U,S,V,temp;
        Index feedbacktimebin = findtype(psi.A(timestep+1-Nfb),Site);
        for(int i=timestep+1-Nfb; i<timestep; i++)
        {
            temp = psi.A(i)*psi.A(i+1);
            U = ITensor(feedbacktimebin, commonIndex(psi.A(i+1),psi.A(i+2)));
            svd(temp,U,S,V,args);
            psi.setA(i,V);
            psi.setA(i+1,U*S);
        }
    }
    
    void TLS_Liouville::swap_feedbackbin_to_origin(int timestep)
    {
        ITensor temp,U,S,V;
        // feedbackbin carries the index with a difference of Nfb to system, thus its index is (timestep+1)-Nfb (the +1 because we start counting at Nfb)
        Index feedbacktimebin = findtype(psi.A(timestep),Site);
        // loop ends when feedback bin is located one bin right of its original position
        // because in the last step, the orthocenter needs to be placed
        // at the future feedback bin for the next time step
        for(int i=timestep; i>timestep-Nfb+2; i--) 
        {
            temp = psi.A(i-1)*psi.A(i);
            U = ITensor(feedbacktimebin, leftLinkInd(psi,i-1));
            svd(temp,U,S,V,args);
            psi.setA(i,V);
            psi.setA(i-1,U*S);
        }
        temp = psi.A(timestep-Nfb+1)*psi.A(timestep-Nfb+2); 
        U = ITensor(feedbacktimebin, leftLinkInd(psi,timestep-Nfb+1));
        svd(temp,U,S,V,args);
        psi.setA(timestep-Nfb+2,V*S);
        psi.setA(timestep-Nfb+1,U);
    }
    
    void TLS_Liouville::timeevolve_system_with_feedback()
    {
        println("\n############################################################## \n \n");
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "a");
        TLS_Liouville::setup_indexbin();
        TLS_Liouville::initialize_MPS_with_feedback();
        TLS_Liouville::set_up_Uevo_with_feedback();
        println("starting timeevolution of system... \n");
        printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", 0.0*dt, t_end*dt, TLS_Liouville::excited_TLS_occupation(psi.A(1), Site),  TLS_Liouville::trace(psi.A(1), Site) ,TLS_Liouville::rhosquare(psi.A(1), Site));
        if (use_feedback=="yes")
        {
            println("Using feedback!");
        }
        else if (use_feedback=="no")
        {
            println("Without feedback!");
        }
        else
        {
            println("Could not identify use of feedback. Exiting");
            throw std::exception();
        }
        for(int j=Nfb; j<ttotal; j++)
        {
            print(psi);
            TLS_Liouville::apply_Uevo_with_feedback(j);
            if(j<(ttotal-1))
            {
                // Uevo needs to act on future bin of MPS for next timestep:
                Uevo *= delta(bin[j+1],bin[j+2]);
                Uevo *= delta(prime(bin[j+1]),prime(bin[j+2]));
                Uevo *= delta(bin[j-Nfb+1],bin[j-Nfb+2]);
                Uevo *= delta(prime(bin[j-Nfb+1]),prime(bin[j-Nfb+2]));
            }
        }
        // system is now at the end of MPS, thus at psi.A(ttotal+1).
        // last bin with feedback interaction has index bin[Nfb] and is at position psi.A(Nfb).
        fclose(datastream_timeevo);
        println("\nfinished timeevolution of system.\n \n");
    }
    
    void TLS_Liouville::calc_correlation()
    {
        println("start calculating g(1) and g(2) functions\n");
        datastream_correlation = fopen(file_identifier_correlation.c_str(), "a");
        // set up reference annihilator:
        // first find reference bin: this is the last bin
        // which has been in feedback-interaction with the system. 
        // after time loop, this bin is located at psi.A(ttotal-Nfb)
        // and the orthocenter is at psi.A(ttotal-Nfb+1)
        // so move orthocenter to reference bin first, but do not move the bins!
        ITensor U, S, V, temp;
        referencetimebin = ttotal-Nfb;
        temp = psi.A(ttotal-Nfb)*psi.A(ttotal-Nfb+1);
        U = ITensor(bin[referencetimebin], leftLinkInd(psi,ttotal-Nfb));
        svd(temp,U,S,V,args);
        psi.setA(ttotal-Nfb,U*S);
        psi.setA(ttotal-Nfb+1,V);
        Sm_RefBin = TLS_Liouville::annihilator(bin[referencetimebin]);
        
        // allocate vectors
        g1_bath = std::vector<Complex>(spectrum_steps+2);
        g2_bath = std::vector<Complex>(spectrum_steps+2);
        pop_bath = std::vector<Complex>(spectrum_steps+2);
        
        // calculate self correlation
        g1_tensor = dag(noprime(psi.A(referencetimebin)*Sm_RefBin)) * noprime(Sm_RefBin*psi.A(referencetimebin));
        g1_bath[0] = g1_tensor.cplx();
        
        g2_tensor = dag(noprime(Sm_RefBin*noprime(Sm_RefBin*psi.A(referencetimebin)))) *
        noprime(Sm_RefBin*noprime(Sm_RefBin*psi.A(referencetimebin)));
        g2_bath[0] = g2_tensor.cplx();

        pop_bath_tensor = dag(noprime(Sm_RefBin*psi.A(referencetimebin))) * noprime(Sm_RefBin*psi.A(referencetimebin));
        pop_bath[0] = pop_bath_tensor.cplx();
    
        for(int i=referencetimebin;i>(referencetimebin-spectrum_steps);i--)
        {
            // calculate g1, g2 and bath population
            int k = referencetimebin-i+1;
            Index s = findtype(psi.A(i-1), Site);
            temp = psi.A(i)*psi.A(i-1);
            g1_bath[k] = TLS_Liouville::calc_g1_bath(s, temp);
            g2_bath[k] = TLS_Liouville::calc_g2_bath(s, temp);
            pop_bath[k] = TLS_Liouville::calc_pop_bath(temp);
            
            if((i%save_nth_timestep)==0) 
            {
                fprintf(datastream_correlation,"%e \t %e \t %e \t %e \t %e \t %e \n",dt*k,g1_bath[k].real(),g1_bath[k].imag(),g2_bath[k].real(), g2_bath[k].imag(), pop_bath[k].real()*pop_bath[k].real());
            }
            
            if((i%print_nth_timestep)==0) 
            {
                printf("g1_bath[%i]=%e+i%e -- g2_bath[%i]=%e+i%e -- (bath_pop)^2=%e \n", k, g1_bath[k].real(), g1_bath[k].imag(), k, g2_bath[k].real(), g2_bath[k].imag(), pop_bath[k].real()*pop_bath[k].real());
            }
            
            // swap bin[referencetimebin] and orthocenter one step to the left for next step of calculation
            if(i>(referencetimebin-spectrum_steps+1))
            {
                ITensor U1, V1, S1;
                U1 = ITensor(bin[referencetimebin],commonIndex(psi.A(i-1), psi.A(i-2)));
                svd(temp,U1,S1,V1);
                psi.setA(i-1,U1*S1);
                psi.setA(i,V1);
            }
        }
        fclose(datastream_correlation);
        println("\nfinished calculating g(1) and g(2) functions\n");
    }
    
    void TLS_Liouville::calc_spectrum()
    {
        println("start calculating fourier trafo for spectrum....");
        // params for normalization and calculation
        spectrum = std::vector<Complex>(spectrum_steps+2);
        dw = spectrum_interval/spectrum_steps;
        
        // needed to calculate the spectrum without coherent part
        std::complex<double> b_offset = g1_bath[spectrum_steps]; 
        // needed for normalization TODO correct would be to use (bath_pop)^2, so check
        double g2_norm = g2_bath[spectrum_steps].real(); 
        
        datastream_spectrum = fopen(file_identifier_spectrum.c_str(), "a");
        
        // calculate spectral density with Wiener-Chinschin
        // thus make a fourier trafo of g(1)
        for(int i=0;i<=spectrum_steps;i++)
        {
            om = -0.5*(spectrum_interval)+i*dw;
            spectrum[i] = 0.0;
            for(int k=0;k<=spectrum_steps;k++)
            {
                spectrum[i] += exp(-Cplx_i*om*k*dt)*(g1_bath[k]-b_offset); 
            }
            fprintf(datastream_spectrum,"%e \t %e \t %e \t %e \n", om, spectrum[i].real(), spectrum[i].imag(), g2_bath[i].real()/g2_norm);
        }
        fclose(datastream_spectrum);
        println("...done\n");
    }
    

}//eof namespace




