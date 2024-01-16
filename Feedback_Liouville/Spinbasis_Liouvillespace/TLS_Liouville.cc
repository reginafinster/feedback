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
        GammaOut = input.getReal("GammaOut");
        GammaIn = input.getReal("GammaIn");
        Kappa = input.getReal("Kappa");
        
        // params for feedback
        Nfb = input.getInt("feedback_time");
        phi = input.getReal("feedbackphase_in_units_of_phi",0.0);
        
        // params for spectrum
        spectrum_steps = input.getInt("spectrum_steps");
        spectrum_interval = input.getReal("spectrum_interval");  
        
        // params for program control
        save_nth_timestep = input.getInt("save_nth_timestep");
        print_nth_timestep = input.getInt("print_nth_timestep");
        
        // calculate some params
        ttotal = t_end + Nfb;
        phi *= 3.141592654;
        
        GammaOut = pow(GammaOut,2.0);
        GammaIn = pow(GammaIn,2.0);
        
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

    double TLS_Liouville::excited_occupation(ITensor psitensor)
    {
        Index s = findtype(psitensor, Site);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(1), prime(s(1)), 0.5);
        SpSm.set(s(4), prime(s(1)), -0.5);
        ITensor temp = (SpSm*psitensor)*(SpSm*psitensor);
        double occ = temp.real();
        occ /= pow(2,ttotal);
        occ = pow(occ,0.5);
        return occ;
    }
    
    double TLS_Liouville::ground_occupation(ITensor psitensor)
    {
        Index s = findtype(psitensor, Site);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(1), prime(s(1)), 0.5);
        SpSm.set(s(4), prime(s(1)), 0.5);
        ITensor temp = (SpSm*psitensor)*(SpSm*psitensor);
        double occ = temp.real();
        occ /= pow(2,ttotal);
        occ = pow(occ,0.5);
        return occ;
    }
    
    double TLS_Liouville::trace(ITensor psitensor)
    {
        Index s = findtype(psitensor, Site);
        ITensor TraceOp = ITensor(s, prime(s));
        TraceOp.set(s(1), prime(s(1)), 1.0);
        ITensor temp = (TraceOp*psitensor) * (TraceOp*psitensor);
        double trace = temp.real();
        trace /= pow(2,ttotal);
        trace = pow(trace,0.5);
        return trace;
    }

    double TLS_Liouville::rhosquare(ITensor psi)
    {
        double rhosq = (dag(psi)*psi).real();
        rhosq /= pow(2,ttotal+1);
        return rhosq;
    }
    
    std::complex<double> TLS_Liouville::calc_g1_bath(Index refbin, Index timebin, ITensor temp)
    {
//         ITensor g1_operator = ITensor(refbin, prime(refbin), timebin, prime(timebin));
//         g1_operator.set(refbin(2), prime(refbin(2)), timebin(3), prime(timebin)(3), 1.0);
//         g1_tensor = dag(prime(temp, Site)) * (g1_operator*temp);
//         std::complex<double> g1_cplx = g1_tensor.cplx();
         std::complex<double> g1_cplx = 0.0;
         return g1_cplx;
    }
    
    std::complex<double> TLS_Liouville::calc_g2_bath(Index refbin, Index timebin, ITensor temp)
    {
        std::complex<double> g1_cplx = 0.0;
        return g1_cplx;
//         // set up sigma- and sigma+ for the subspace of the timebin
//         ITensor Sm_TimeBin = ITensor(timebin, prime(timebin));
//         Sm_TimeBin.set(timebin(3), prime(timebin)(1), 1.0);
//         Sm_TimeBin.set(timebin(4), prime(timebin)(2), 1.0);
//         
//         ITensor Sp_TimeBin = ITensor(timebin, prime(timebin));
//         Sp_TimeBin.set(timebin(1), prime(timebin)(3), 1.0);
//         Sp_TimeBin.set(timebin(2), prime(timebin)(4), 1.0);
//         
//         // calculate g(2) in density matrix picture:
//         ITensor g2_tensor;
//         g2_tensor = noprime(Sp_RefBin*noprime(Sp_TimeBin*(noprime(Sm_TimeBin*(noprime(Sm_RefBin*temp))))));
        
//         Print(temp);
        
//         // set up operator:
//         Sm_RefBin = ITensor(bin[referencetimebin], prime(bin[referencetimebin]));
//         ITensor g2_operator = ITensor(refbin, prime(refbin), bin[referencetimebin-2], prime(bin[referencetimebin-2]), timebin, prime(timebin));
//         g2_operator.set(refbin(4), prime(refbin(4)), timebin(4), bin[referencetimebin-2](4), prime(bin[referencetimebin-2])(4), prime(timebin)(4), 1.0);
//       
//         // calculate g(2) in density matrix picture:
//         ITensor g2_tensor;
//         PrintData(g2_operator);
//         g2_tensor = dag(temp)*noprime(g2_operator*temp);
//         double g2 = g2_tensor.real();
//         g2 = pow(g2,0.5);
        
        
//         // calculate g(2) in density matrix picture:
//         ITensor g2_tensor;
//         g2_tensor = noprime(g2_operator*temp);
//         PrintData(g2_tensor);
        
//         ITensor TraceOp1 = ITensor(refbin, prime(refbin), timebin, prime(timebin));
//         TraceOp1.set(refbin(1), prime(refbin)(1), timebin(1), prime(timebin)(1), 1.0);
//         ITensor temp1 = dag(prime(g2_tensor, Site)) * (TraceOp1*g2_tensor);
//         double trace1 = temp1.real();
//         trace1 = pow(trace1,0.5);
//         
//         ITensor TraceOp2 = ITensor(refbin, prime(refbin), timebin, prime(timebin));
//         TraceOp2.set(refbin(1), prime(refbin)(1), timebin(4), prime(timebin)(4), 1.0);
//         ITensor temp2 = dag(prime(g2_tensor, Site)) * (TraceOp2*g2_tensor);
//         double trace2 = temp2.real();
//         trace2 = pow(trace2,0.5);
//         
//         ITensor TraceOp3 = ITensor(refbin, prime(refbin), timebin, prime(timebin));
//         TraceOp3.set(refbin(4), prime(refbin)(4), timebin(1), prime(timebin)(1), 1.0);
//         ITensor temp3 = dag(prime(g2_tensor, Site)) * (TraceOp3*g2_tensor);
//         double trace3 = temp3.real();
//         trace3 = pow(trace3,0.5);
//         
//         ITensor TraceOp4 = ITensor(refbin, prime(refbin), timebin, prime(timebin));
//         TraceOp4.set(refbin(4), prime(refbin)(4), timebin(4), prime(timebin)(4), 1.0);
//         ITensor temp4 = dag(prime(g2_tensor, Site)) * (TraceOp4*g2_tensor);
//         Print(temp4.cplx());
//         double trace4 = temp4.real();
//         trace4 = pow(trace4,0.5);
        
//         Print(trace1);
//         Print(trace2);
//         Print(trace3);
//         Print(trace4);
//         double g2 = trace1 + trace2 + trace3 + trace4;
//        
//         return g2;
    }
    
    void TLS_Liouville::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::string date = TLS_Liouville::get_date_as_string();
        std::ostringstream oss1, oss2, oss3;
        oss1 << "out/timeevo_" << date.c_str() << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamOut"<< GammaOut << "_GamIn" << GammaIn << "_Kappa" << Kappa << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phi" << (phi/3.141592654) << ".dat";
        file_identifier_timeevo = oss1.str();
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "w");
        fprintf(datastream_timeevo,"#time \tS+S- \tS-S+ \t S+S- + S-S+ \t trace \t purity \n");
        fclose(datastream_timeevo);
        
        oss2 << "out/corr_" << date.c_str() << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamOut"<< GammaOut << "_GamIn" << GammaIn << "_Kappa" << Kappa << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phi" << (phi/3.141592654) << ".dat";
        file_identifier_correlation = oss2.str();
        datastream_correlation = fopen(file_identifier_correlation.c_str(), "w");
        fprintf(datastream_correlation, "#time \t Re[g1] \t Im[g1] \t Re[g2] \t Im[g2] \t pop_bath \t pop_bath^2 \n");
        fclose(datastream_correlation);
        
        oss3 << "out/spec_" << date.c_str() << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamOut"<< GammaOut << "_GamIn" << GammaIn << "_Kappa" << Kappa << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phi" << (phi/3.141592654) << ".dat";
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
            MPSTensor.set(bin[i](4),1.0);
            //MPSTensor *= (0.5*pow(2,0.5));
            psi.setA((i),MPSTensor);
        }
        
        MPSTensor = ITensor(bin[0]);
        if(init_cond=="ground")
        {
            MPSTensor.set(bin[0](1),1.0);
            MPSTensor.set(bin[0](4),1.0);
            //MPSTensor *= (0.5*pow(2,0.5));
            println("Initialize system in ground state\n");
        }
        else if(init_cond=="excited")
        {
            MPSTensor.set(bin[0](1),1.0);
            MPSTensor.set(bin[0](4),-1.0);
            //MPSTensor *= (0.5*pow(2,0.5));
            println("Initialize system in excited state\n");
        }
        psi.setA((Nfb+1),MPSTensor);
        
        for (int i=Nfb+2;i<ttotal+2;i++)
        {
            MPSTensor = ITensor(bin[i-1]);
            MPSTensor.set(bin[i-1](1),1.0);
            MPSTensor.set(bin[i-1](4),1.0);
            //MPSTensor *= (0.5*pow(2,0.5));
            psi.setA((i),MPSTensor);
        }
        // orthogonalize in order to establish links
        // orthocenter is now at psi.A(1)
        psi.orthogonalize();
    }
    
    void TLS_Liouville::set_up_Uevo_with_feedback()
    {
        Lsys = ITensor(bin[0],prime(bin[0]),bin[Nfb+1],prime(bin[Nfb+1]));
        for (int p=1; p<5; p++)
        {
            Lsys.set(bin[0](3),bin[Nfb+1](p),prime(bin[0](4)),prime(bin[Nfb+1](p)),2.0*Omega*dt);
            Lsys.set(bin[0](4),bin[Nfb+1](p),prime(bin[0](3)),prime(bin[Nfb+1](p)),-2.0*Omega*dt);
        }
        
        Lint = ITensor(bin[0],prime(bin[0]),bin[Nfb+1],prime(bin[Nfb+1]));
        
        Lint.set(bin[0](1),bin[Nfb+1](2),prime(bin[0](2)),prime(bin[Nfb+1](4)), -Gamma*pow(dt,0.5));
        Lint.set(bin[0](1),bin[Nfb+1](3),prime(bin[0](3)),prime(bin[Nfb+1](4)), -Gamma*pow(dt,0.5));
        Lint.set(bin[0](1),bin[Nfb+1](4),prime(bin[0](2)),prime(bin[Nfb+1](2)), +Gamma*pow(dt,0.5));
        Lint.set(bin[0](1),bin[Nfb+1](4),prime(bin[0](3)),prime(bin[Nfb+1](3)), +Gamma*pow(dt,0.5));

        Lint.set(bin[0](2),bin[Nfb+1](1),prime(bin[0](4)),prime(bin[Nfb+1](2)), +Gamma*pow(dt,0.5));
        Lint.set(bin[0](2),bin[Nfb+1](2),prime(bin[0](4)),prime(bin[Nfb+1](1)), +Gamma*pow(dt,0.5));
        Lint.set(bin[0](2),bin[Nfb+1](2),prime(bin[0](1)),prime(bin[Nfb+1](4)), -Gamma*pow(dt,0.5));
        Lint.set(bin[0](2),bin[Nfb+1](4),prime(bin[0](1)),prime(bin[Nfb+1](2)), +Gamma*pow(dt,0.5));

        Lint.set(bin[0](3),bin[Nfb+1](1),prime(bin[0](4)),prime(bin[Nfb+1](3)), +Gamma*pow(dt,0.5));
        Lint.set(bin[0](3),bin[Nfb+1](3),prime(bin[0](4)),prime(bin[Nfb+1](1)), +Gamma*pow(dt,0.5));
        Lint.set(bin[0](3),bin[Nfb+1](3),prime(bin[0](1)),prime(bin[Nfb+1](4)), -Gamma*pow(dt,0.5));
        Lint.set(bin[0](3),bin[Nfb+1](4),prime(bin[0](1)),prime(bin[Nfb+1](3)), +Gamma*pow(dt,0.5));

        Lint.set(bin[0](4),bin[Nfb+1](1),prime(bin[0](3)),prime(bin[Nfb+1](3)), -Gamma*pow(dt,0.5));
        Lint.set(bin[0](4),bin[Nfb+1](1),prime(bin[0](2)),prime(bin[Nfb+1](2)), -Gamma*pow(dt,0.5));
        Lint.set(bin[0](4),bin[Nfb+1](2),prime(bin[0](2)),prime(bin[Nfb+1](1)), -Gamma*pow(dt,0.5));
        Lint.set(bin[0](4),bin[Nfb+1](3),prime(bin[0](3)),prime(bin[Nfb+1](1)), -Gamma*pow(dt,0.5));
        
        
        Lfb = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));

        Lfb.set(bin[0](1),bin[1](2),prime(bin[0](2)),prime(bin[1](4)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](1),bin[1](3),prime(bin[0](3)),prime(bin[1](4)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](1),bin[1](4),prime(bin[0](2)),prime(bin[1](2)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](1),bin[1](4),prime(bin[0](3)),prime(bin[1](3)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));

        Lfb.set(bin[0](2),bin[1](1),prime(bin[0](4)),prime(bin[1](2)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](2),bin[1](2),prime(bin[0](4)),prime(bin[1](1)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](2),bin[1](2),prime(bin[0](1)),prime(bin[1](4)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](2),bin[1](4),prime(bin[0](1)),prime(bin[1](2)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));

        Lfb.set(bin[0](3),bin[1](1),prime(bin[0](4)),prime(bin[1](3)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](3),bin[1](3),prime(bin[0](4)),prime(bin[1](1)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](3),bin[1](3),prime(bin[0](1)),prime(bin[1](4)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](3),bin[1](4),prime(bin[0](1)),prime(bin[1](3)), +Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));

        Lfb.set(bin[0](4),bin[1](1),prime(bin[0](3)),prime(bin[1](3)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](4),bin[1](1),prime(bin[0](2)),prime(bin[1](2)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](4),bin[1](2),prime(bin[0](2)),prime(bin[1](1)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        Lfb.set(bin[0](4),bin[1](3),prime(bin[0](3)),prime(bin[1](1)), -Gamma*0.5*(exp(Cplx_i*phi)+exp(-Cplx_i*phi))*pow(dt,0.5));
        
//         LindbladLoss = ITensor(bin[0],prime(bin[0]));
//         LindbladLoss.set(bin[0](4),prime(bin[0](1)), 2.0*GammaOut*dt);
//         LindbladLoss.set(bin[0](2),prime(bin[0](2)), -1.0*GammaOut*dt);
//         LindbladLoss.set(bin[0](3),prime(bin[0](3)), -1.0*GammaOut*dt);
//         LindbladLoss.set(bin[0](4),prime(bin[0](4)), -2.0*GammaOut*dt);
//         
//         LindbladPump = ITensor(bin[0],prime(bin[0]));
//         LindbladPump.set(bin[0](1),prime(bin[0](4)), 2.0*GammaIn*dt);
//         LindbladPump.set(bin[0](1),prime(bin[0](1)), -2.0*GammaIn*dt);
//         LindbladPump.set(bin[0](2),prime(bin[0](2)), -1.0*GammaIn*dt);
//         LindbladPump.set(bin[0](3),prime(bin[0](3)), -1.0*GammaIn*dt);
//         
//         PureDephasing = ITensor(bin[0],prime(bin[0]));
//         PureDephasing.set(bin[0](2), prime(bin[0](2)), -1.0*Kappa*dt);
//         PureDephasing.set(bin[0](3), prime(bin[0](3)), -1.0*Kappa*dt);
/*        
        Lindblad = ITensor(bin[0],prime(bin[0]));
        Lindblad = LindbladLoss + LindbladPump + PureDephasing;*/
        
        ITensor identity_system = TLS_Liouville::identity(bin[0]);
        ITensor identity_int = TLS_Liouville::identity(bin[Nfb+1]);
        ITensor identity_fb = TLS_Liouville::identity(bin[1]);
        
//         Lindblad *= (identity_int*identity_fb);
        Lsys *= identity_fb;
        Lfb *= identity_int;
        Lint *= identity_fb;
        
//         ITensor U_fb_1  = Lsys + Lint + Lfb + Lindblad; 
        ITensor U_fb_1  = Lsys + Lint + Lfb; 
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
        fprintf(datastream_timeevo, "%f \t %e \t %e \t %e \t %e \t %e \n", (timestep-Nfb)*dt, TLS_Liouville::excited_occupation(psi.A(timestep+1)), TLS_Liouville::ground_occupation(psi.A(timestep+1)), (TLS_Liouville::excited_occupation(psi.A(timestep+1))+TLS_Liouville::ground_occupation(psi.A(timestep+1))), TLS_Liouville::trace(psi.A(timestep+1)), TLS_Liouville::rhosquare(psi.A(timestep+1)));
        if ((timestep%10)==0) 
        {
            printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", (timestep-Nfb)*dt, t_end*dt, TLS_Liouville::excited_occupation(psi.A(timestep+1)),  TLS_Liouville::trace(psi.A(timestep+1)) ,TLS_Liouville::rhosquare(psi.A(timestep+1)));
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
        printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", 0.0*dt, t_end*dt, TLS_Liouville::excited_occupation(psi.A(1)),  TLS_Liouville::trace(psi.A(1)) ,TLS_Liouville::rhosquare(psi.A(1)));
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
        // and the orthocenter is at psi.A(ttotal-Nfb+1) // TODO das hab ich überprüft, das stimmt!
        // so move orthocenter to reference bin first, but do not move the bins!
        ITensor U, S, V, temp;
        // if spectrum of feedback bath is required: referencetimebin = ttotal-Nfb.
        // if spectrum of interaction without feedback bath: referencetimebin = ttotal.
        // system is at the end of MPS after time evolution, thus at psi.A(ttotal+1)
        // thus last bath bin is at psi.A(ttotal).
        referencetimebin = ttotal-Nfb;
        temp = psi.A(ttotal-Nfb)*psi.A(ttotal-Nfb+1);
        U = ITensor(bin[referencetimebin], leftLinkInd(psi,ttotal-Nfb));
        svd(temp,U,S,V,args);
        psi.setA(ttotal-Nfb,U*S);
        psi.setA(ttotal-Nfb+1,V);
        
        // set up sigma- and sigma+ for the reference timebin
        Sm_RefBin = ITensor(bin[referencetimebin], prime(bin[referencetimebin]));
        Sm_RefBin.set(bin[referencetimebin](3), prime(bin[referencetimebin])(1), 1.0);
        Sm_RefBin.set(bin[referencetimebin](4), prime(bin[referencetimebin])(2), 1.0);
        
        Sp_RefBin = ITensor(bin[referencetimebin], prime(bin[referencetimebin]));
        Sp_RefBin.set(bin[referencetimebin](1), prime(bin[referencetimebin])(3), 1.0);
        Sp_RefBin.set(bin[referencetimebin](2), prime(bin[referencetimebin])(4), 1.0);
        
        // allocate vectors
        g1_bath = std::vector<Complex>(spectrum_steps+2);
        g2_bath = std::vector<Complex>(spectrum_steps+2);
        pop_bath = std::vector<double>(spectrum_steps+2);
        
        // calculate self correlation
        g1_bath[0] = TLS_Liouville::excited_occupation(psi.A(referencetimebin));
        g2_bath[0] = 0.0;
        pop_bath[0] = TLS_Liouville::excited_occupation(psi.A(referencetimebin));
        
        // contract rho
        temp = psi.A(1);
        for (int j = 2; j<=t_end; j++)
        {
            temp *= psi.A(j);
        }
    
        for(int i=referencetimebin;i>(referencetimebin-spectrum_steps);i--)
        {
            // calculate g1, g2 and bath population
            int k = referencetimebin-i+1;
            Index s = findtype(psi.A(i-1), Site);
            //temp = psi.A(i)*psi.A(i-1);
        
            //g1_bath[k] = TLS_Liouville::calc_g1_bath(bin[referencetimebin], s, temp);
            g2_bath[k] = TLS_Liouville::calc_g2_bath(bin[referencetimebin], s, temp);
            pop_bath[k] = TLS_Liouville::excited_occupation(psi.A(i));
            
            if((i%save_nth_timestep)==0) 
            {
                fprintf(datastream_correlation,"%e \t %e \t %e \t %e \t %e \t %e \t %e \n", dt*k, g1_bath[k].real(), g1_bath[k].imag(), g2_bath[k].real() , g2_bath[k].imag(), pop_bath[k], pop_bath[k]*pop_bath[k]);
            }
            
            if((i%print_nth_timestep)==0) 
            {
                printf("g1_bath[%i]=%e+i%e -- g2_bath[%i]=%e+i%e -- (bath_pop)^2=%e \n", k, g1_bath[k].real(), g1_bath[k].imag(), k, g2_bath[k].real(), g2_bath[k].imag(), pop_bath[k]*pop_bath[k]);
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




