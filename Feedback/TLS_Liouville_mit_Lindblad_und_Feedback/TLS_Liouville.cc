#include "TLS_Liouville.h"


namespace itensor{
    
    TLS_Liouville::TLS_Liouville(std::string infile)
    {
        InputGroup input = InputGroup(infile,"input");
        use_feedback = input.getString("use_feedback");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        dt = input.getReal("time_step");
        t_end = input.getReal("time_end");
        Omega = input.getReal("Omega_TLS");
        Gamma = input.getReal("Gamma");
        LindbladGamma = input.getReal("LindbladGamma");
        LindbladGamma /= 3.141592654;
        LindbladGamma = pow(LindbladGamma,2.0);
        LindbladGamma *= 2.0;
        Nfb = input.getInt("feedback_time");
        init_cond = input.getString("initial_condition");
        args = Args("Cutoff=",cutoff,"Maxm=",maxm);
        ttotal = t_end/dt;
        phi = input.getReal("feedbackphase_in_units_of_phi",0.0);
        phi *= 3.141592654;
        if (use_feedback=="yes")
        {
            ttotal += Nfb;
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

    double TLS_Liouville::rhosquare(ITensor psi, IndexType type)
    {
        return (dag(psi)*psi).real();
    }

    void TLS_Liouville::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::string date = TLS_Liouville::get_date_as_string();
        std::ostringstream oss;
        oss << "out/TLS_" << date.c_str() << "feedback_" << use_feedback << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamLind"<< LindbladGamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier = oss.str();
        datastream = fopen(file_identifier.c_str(), "w");
        fprintf(datastream,"#time \tS+S- \tS-S+ \t trace \t purity \n");
        fclose(datastream);
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
    
    void TLS_Liouville::initialize_MPS_without_feedback()
    {
        // We need one tensor for each timestep and first one is system
        // system may be in ground or exited state, timebins are initialized in vacuum
        psi = MPS(ttotal+1);
        ITensor MPSTensor = ITensor();
        
        MPSTensor = ITensor(bin[0]);
        if(init_cond=="ground")
        {
            MPSTensor.set(bin[0](1),1.0);
            println("Initialize system in ground state\n");
        }
        else if(init_cond=="excited")
        {
            MPSTensor.set(bin[0](4),1.0);
            println("Initialize system in excited state\n");
        }
        psi.setA((1),MPSTensor);

        for (int i=1;i<ttotal+1;i++)
        {
            MPSTensor = ITensor(bin[i]);
            MPSTensor.set(bin[i](1),1.0);
            psi.setA((i+1),MPSTensor);
        }
        
        psi.orthogonalize();
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
        psi.orthogonalize();
    }
    
    void TLS_Liouville::set_up_Uevo_without_feedback()
    {
        
        Lsys = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));
        for (int p=1; p<5; p++)
        {
            Lsys.set(bin[0](1),bin[1](p),prime(bin[0](3)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](2),bin[1](p),prime(bin[0](4)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](3),bin[1](p),prime(bin[0](1)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](4),bin[1](p),prime(bin[0](2)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](2),bin[1](p),prime(bin[0](1)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](4),bin[1](p),prime(bin[0](3)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](1),bin[1](p),prime(bin[0](2)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](3),bin[1](p),prime(bin[0](4)),prime(bin[1](p)),-Cplx_i*Omega*dt);
        }
        
        Lint = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));
        
        for (int p=1; p<=2; p++)
        {
            Lint.set(bin[0](1),bin[1](p+2),prime(bin[0](3)),prime(bin[1](p)), Gamma*pow(dt,0.5));
            Lint.set(bin[0](2),bin[1](p+2),prime(bin[0](4)),prime(bin[1](p)), Gamma*pow(dt,0.5));
        }
        
        for (int p=3; p<=4; p++)
        {
            Lint.set(bin[0](3),bin[1](p-2),prime(bin[0](1)),prime(bin[1](p)), -Gamma*pow(dt,0.5));
            Lint.set(bin[0](4),bin[1](p-2),prime(bin[0](2)),prime(bin[1](p)), -Gamma*pow(dt,0.5));
        }
        
        for (int p=2; p<=4; p+=2)
        {
            Lint.set(bin[0](2),bin[1](p-1),prime(bin[0](1)),prime(bin[1](p)), -Gamma*pow(dt,0.5));
            Lint.set(bin[0](4),bin[1](p-1),prime(bin[0](3)),prime(bin[1](p)), -Gamma*pow(dt,0.5));
        }
        
        for (int p=1; p<=3; p+=2)
        {
            Lint.set(bin[0](1),bin[1](p+1),prime(bin[0](2)),prime(bin[1](p)), Gamma*pow(dt,0.5));
            Lint.set(bin[0](3),bin[1](p+1),prime(bin[0](4)),prime(bin[1](p)), Gamma*pow(dt,0.5));
        }
        
        // Levo = identity_int*identity_system + Lsys + Lint + 0.5*mapprime(prime(Lint)*Lint,2,1);

        Lindblad = ITensor(bin[0],prime(bin[0]));
        Lindblad.set(bin[0](4),prime(bin[0](1)), LindbladGamma);
        Lindblad.set(bin[0](2),prime(bin[0](2)), -0.5*LindbladGamma);
        Lindblad.set(bin[0](3),prime(bin[0](3)), -0.5*LindbladGamma);
        Lindblad.set(bin[0](4),prime(bin[0](4)), -1.0*LindbladGamma);
        
        ITensor identity_system = TLS_Liouville::identity(bin[0]);
        ITensor identity_int = TLS_Liouville::identity(bin[1]);
        
        Lindblad *= identity_int;
        
        ITensor U_fb_1  = Lsys + Lint + Lindblad; 
        ITensor U_fb_2  = (1./2.)  * mapprime(U_fb_1*prime(U_fb_1),2,1);
        ITensor U_fb_3  = (1./3.)  * mapprime(U_fb_1*prime(U_fb_2),2,1);
        ITensor U_fb_4  = (1./4.)  * mapprime(U_fb_1*prime(U_fb_3),2,1);
        ITensor U_fb_5  = (1./5.)  * mapprime(U_fb_1*prime(U_fb_4),2,1);
        ITensor U_fb_6  = (1./6.)  * mapprime(U_fb_1*prime(U_fb_5),2,1);
        ITensor U_fb_7  = (1./7.)  * mapprime(U_fb_1*prime(U_fb_6),2,1);
        ITensor U_fb_8  = (1./8.)  * mapprime(U_fb_1*prime(U_fb_7),2,1);
        ITensor U_fb_9  = (1./9.)  * mapprime(U_fb_1*prime(U_fb_8),2,1);
        ITensor U_fb_10 = (1./10.) * mapprime(U_fb_1*prime(U_fb_9),2,1);
    
        Uevo = identity_system*identity_int + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7 + U_fb_8 + U_fb_9 + U_fb_10;
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
                
        Lindblad = ITensor(bin[0],prime(bin[0]));
        Lindblad.set(bin[0](4),prime(bin[0](1)), LindbladGamma);
        Lindblad.set(bin[0](2),prime(bin[0](2)), -0.5*LindbladGamma);
        Lindblad.set(bin[0](3),prime(bin[0](3)), -0.5*LindbladGamma);
        Lindblad.set(bin[0](4),prime(bin[0](4)), -1.0*LindbladGamma);
        
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
    
    void TLS_Liouville::apply_Uevo_without_feedback(int timestep)
    {
        ITensor U,S,V,temp;    
        temp = noprime(Uevo*psi.A(timestep+1)*psi.A(timestep+2));
        U = ITensor(bin[timestep+1], leftLinkInd(psi,timestep+1) );
        V = ITensor(bin[0],rightLinkInd(psi,timestep+2));
        svd(temp,U,S,V);
        psi.setA(timestep+1,U);
        psi.setA(timestep+2,V*S);
    }
    
    void TLS_Liouville::apply_Uevo_with_feedback(int timestep)
    {
        TLS_Liouville::swap_feedbackbin_to_system(timestep);
        
        ITensor U,S,V,temp;
        //shift ortho-center from feedbackbin to systembin to compute observable
        temp = psi.A(timestep)*psi.A(timestep+1);
        U = ITensor(bin[0],rightLinkInd(psi,timestep+1));
        svd(temp,U,S,V,args);
        psi.setA(timestep,V);
        psi.setA(timestep+1,U*S);
        fprintf(datastream, "%f \t %e \t %e \t %e \t %e \n", (timestep-Nfb)*dt, TLS_Liouville::excited_TLS_occupation(psi.A(timestep+1), Site), TLS_Liouville::ground_TLS_occupation(psi.A(timestep+1), Site), TLS_Liouville::trace(psi.A(timestep+1), Site), TLS_Liouville::rhosquare(psi.A(timestep+1), Site));
        if ((timestep%10)==0) 
        {
            printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", (timestep-Nfb)*dt, t_end, TLS_Liouville::excited_TLS_occupation(psi.A(timestep+1), Site),  TLS_Liouville::trace(psi.A(timestep+1), Site) ,TLS_Liouville::rhosquare(psi.A(timestep+1), Site));
        }
        
        Index currenttimebin = findtype(psi.A(timestep+2),Site);
        
        temp = noprime(psi.A(timestep)*psi.A(timestep+1)*psi.A(timestep+2)*Uevo);
        
        // bring system to A(timestep+2), thus to the left
        U = ITensor(bin[0],leftLinkInd(psi,timestep+3));
        svd(temp,U,S,V,args);
        V *= S;
        psi.setA(timestep+2,U);
            
        // bring current timebin one position to the left (which is middle position)
        // and feedbacktimebin back to its position at A(timestep) 
        U = ITensor(currenttimebin,commonIndex(U,V));
        temp=V;
        svd(temp,U,S,V,args);
        psi.setA(timestep+1,U);
        psi.setA(timestep,V*S);
        
        TLS_Liouville::swap_feedbackbin_to_origin(timestep);
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
    
    void TLS_Liouville::timeevolve_system_without_feedback()
    // past=j as MPS starts at 1 and bin[] at 0
    // system is always at psi(j+1), future at psi(j+2)
    // Uevo acts on system and future timebin
    // thus we contract these two tensors, apply Uevo and decompose
    // with U(futurebin, link to past), V(systembin, link to future)
    // and set orthocenter to V 
    // U becomes new past bin in MPS
    // V becomes new systembin at next timestep in MPS
    {
        datastream = fopen(file_identifier.c_str(), "a");
        TLS_Liouville::setup_indexbin();
        TLS_Liouville::initialize_MPS_without_feedback();
        printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", 0.0*dt, t_end, TLS_Liouville::excited_TLS_occupation(psi.A(1), Site),  TLS_Liouville::trace(psi.A(1), Site) ,TLS_Liouville::rhosquare(psi.A(1), Site));
        
        TLS_Liouville::set_up_Uevo_without_feedback();
        for(int j = 0; j<ttotal; j++)
        {
            fprintf(datastream, "%f \t %e \t %e \t %e \t %e \n", j*dt, TLS_Liouville::excited_TLS_occupation(psi.A(j+1), Site), TLS_Liouville::ground_TLS_occupation(psi.A(j+1), Site), TLS_Liouville::trace(psi.A(j+1), Site), TLS_Liouville::rhosquare(psi.A(j+1), Site));
            if ((j%10)==0) 
            {
                printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", j*dt, t_end, TLS_Liouville::excited_TLS_occupation(psi.A(j+1), Site),  TLS_Liouville::trace(psi.A(j+1), Site) ,TLS_Liouville::rhosquare(psi.A(j+1), Site));
            }
            TLS_Liouville::apply_Uevo_without_feedback(j);
            if(j<(ttotal-1))
            {
                Uevo = Uevo*delta(bin[j+1],bin[j+2])*delta(prime(bin[j+1]),prime(bin[j+2]));
            }
        }
        
        fclose(datastream);
    }
    
    void TLS_Liouville::timeevolve_system_with_feedback()
    {
        datastream = fopen(file_identifier.c_str(), "a");
        TLS_Liouville::setup_indexbin();
        TLS_Liouville::initialize_MPS_with_feedback();
        printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", 0.0*dt, t_end, TLS_Liouville::excited_TLS_occupation(psi.A(1), Site),  TLS_Liouville::trace(psi.A(1), Site) ,TLS_Liouville::rhosquare(psi.A(1), Site));
        TLS_Liouville::set_up_Uevo_with_feedback();
        for(int j=Nfb; j<ttotal; j++)
        {
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
        fclose(datastream);
    }

}//eof namespace




