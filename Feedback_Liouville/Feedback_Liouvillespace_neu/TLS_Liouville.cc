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
        
//        Gamma = pow(Gamma,0.5);
        
//        GammaOut = pow(GammaOut,2.0);
//        GammaIn = pow(GammaIn,2.0);
        
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

    double TLS_Liouville::excited_occupation(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(4), prime(s(4)), 1.0);
        ITensor temp = dag(prime(psitensor, Site)) * (SpSm*psitensor);
        double occ = temp.real();
        occ = pow(occ,0.5);
        return occ;
    }
    
//     double TLS_Liouville::excited_occupation_test(ITensor psitensor, IndexType type)
//     {
//         Index s = findtype(psitensor, type);
//         ITensor SpSm = ITensor(s, prime(s));
//         SpSm.set(s(4), prime(s(4)), 1.0);
//         ITensor temp = SpSm*psitensor;
//         Print(temp);
//         IndexSet inds = temp.inds();
//         PrintData(inds);
//         Print(inds[1].m());
//         ITensor temp_id = ITensor(inds[0], inds[1], inds[2]);
//         for(int i=1; i<=inds[2].m(); i++)
//         {
//             temp_id.set(inds[0](i), inds[1](i), inds[2](i), 1.0);
//         }
//         PrintData(temp_id);
//         temp *= temp_id;
//         PrintData(temp);
//         double occ = temp.real();
//         return occ;
//     }
    
    double TLS_Liouville::ground_occupation(ITensor psitensor, IndexType type)
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
    
    double TLS_Liouville::trace_contracted(ITensor psi_contracted)
    {
        Complex trace_contr = 
        psi_contracted.cplx(bin[0](1),bin[1](1),bin[2](1),bin[3](1)) +
        psi_contracted.cplx(bin[0](4),bin[1](4),bin[2](4),bin[3](4)) +
        
        psi_contracted.cplx(bin[0](4),bin[1](1),bin[2](1),bin[3](1)) +
        psi_contracted.cplx(bin[0](1),bin[1](4),bin[2](1),bin[3](1)) +
        psi_contracted.cplx(bin[0](1),bin[1](1),bin[2](4),bin[3](1)) +
        psi_contracted.cplx(bin[0](1),bin[1](1),bin[2](1),bin[3](4)) +
        
        psi_contracted.cplx(bin[0](4),bin[1](4),bin[2](1),bin[3](1)) +
        psi_contracted.cplx(bin[0](4),bin[1](1),bin[2](4),bin[3](1)) +
        psi_contracted.cplx(bin[0](4),bin[1](1),bin[2](1),bin[3](4)) +
        psi_contracted.cplx(bin[0](1),bin[1](4),bin[2](4),bin[3](1)) +
        psi_contracted.cplx(bin[0](1),bin[1](4),bin[2](1),bin[3](4)) +
        psi_contracted.cplx(bin[0](1),bin[1](1),bin[2](4),bin[3](4)) +
        
        psi_contracted.cplx(bin[0](4),bin[1](4),bin[2](4),bin[3](1)) +
        psi_contracted.cplx(bin[0](1),bin[1](4),bin[2](4),bin[3](4)) +
        psi_contracted.cplx(bin[0](4),bin[1](4),bin[2](1),bin[3](4)) +
        psi_contracted.cplx(bin[0](4),bin[1](1),bin[2](4),bin[3](4));

        return trace_contr.real();
    }
    
    double TLS_Liouville::tr_ex_occ(ITensor psi_contracted, Index s)
    {
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(3), prime(s(3)), 1.0);
        SpSm.set(s(4), prime(s(4)), 1.0);
        ITensor temp = noprime(SpSm*(psi_contracted));
        double ex_occ = TLS_Liouville::trace_contracted(temp);
        return ex_occ;
    }
    
    double TLS_Liouville::tr_gr_occ(ITensor psi_contracted, Index s)
    {
        ITensor SmSp = ITensor(s, prime(s));
        SmSp.set(s(1), prime(s(1)), 1.0);
        SmSp.set(s(2), prime(s(2)), 1.0);
        ITensor temp = noprime(SmSp*(psi_contracted));
        double gr_occ = TLS_Liouville::trace_contracted(temp);
        return gr_occ;
    }
    

    ITensor TLS_Liouville::apply_MPO_for_readout(MPO MPO_Operator)
    {
        int N = ttotal+1;
        // prime site indices of MPS in order to achieve correct multiplication with non-hermitian MPO
        psi.mapprime(0,1,Site);
        
        // declare indices and tensors
        MPS res = MPS(N);
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
            Index si_left = mapprime(commonIndex(psi.A(i), MPO_Operator.A(i), Site),1,0);
            Index si_right = mapprime(commonIndex(psi.A(i+1), MPO_Operator.A(i+1), Site),1,0);

            // contract tensors depending on their position in the chain
            if (i == 1)
            {
                // define additional link indices
                li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
                li_L = commonIndex(MPO_Operator.A(i+1), MPO_Operator.A(i+2), Link);
                // contract Tensors and decompose them again
                *temp = psi.A(i) * MPO_Operator.A(i) * psi.A(i+1) * MPO_Operator.A(i+1);
                *U = ITensor(si_left);
                *V = ITensor(si_right, li_psi, li_L);
            }
            else if (i == (N-1))
            {
                // define additional link indices
                li_res = commonIndex(res.A(i-1), res.A(i), Link);
                // contract Tensors and decompose them again
                *temp = res.A(i) * psi.A(i+1) * MPO_Operator.A(i+1);
                *U = ITensor(si_left, li_res);
                *V = ITensor(si_right);
            }
            else
            {
                // define additional link indices
                li_res = commonIndex(res.A(i-1), res.A(i), Link);
                li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
                li_L = commonIndex(MPO_Operator.A(i+1), MPO_Operator.A(i+2), Link);
                // contract Tensors and decompose them again
                *temp = res.A(i) * psi.A(i+1) * MPO_Operator.A(i+1);
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
        
        // contract rho_MPS
        ITensor contract = res.A(1);
        for (int k = 1; k<N; k++)
        {
            contract *= psi.A(k);
            contract *= res.A(k+1);
        }
        contract *= psi.A(N);
        
        return contract;
    }
    
    double TLS_Liouville::full_trace(ITensor psitensor)
    {
        MPO MPOTrace = MPO(ttotal+1);
        double trace = 0.0;
        ITensor temp;
        
        // dieser Vorgang ist ein Eintrag.
        /////////////////////////////////
        ITensor MPOTensor;
        
        for (int i=1; i<=ttotal+1; i++)
        {
            Index s = findtype(psi.A(i), Site);
            MPOTensor = ITensor(s, prime(s));
            MPOTensor.set(s(1), prime(s(1)),1.0);
            MPOTrace.setA(i, MPOTensor);
        }
        
        MPOTrace.orthogonalize();
        temp = TLS_Liouville::apply_MPO_for_readout(MPOTrace);
        trace += pow(temp.real(), 0.5);
        ///////////////////////////////
        
        return trace;
    }
    
    double TLS_Liouville::full_exc_occ(ITensor psitensor)
    {
        Index psi_s = findtype(psitensor, Site);
        MPO MPOOcc = MPO(ttotal+1);
        ITensor MPOTensor;
        
        for (int i=1; i<=ttotal+1; i++)
        {
            Index s = findtype(psi.A(i), Site);
            if (s==psi_s)
            {
                MPOTensor = ITensor(s, prime(s));
                MPOTensor.set(s(1), prime(s(1)), 0.5);
                MPOTensor.set(s(4), prime(s(1)), -0.5);
                MPOOcc.setA(i, MPOTensor);
            }
            else
            {
                MPOTensor = ITensor(s, prime(s));
                MPOTensor.set(s(1), prime(s(1)),1.0);
                MPOOcc.setA(i, MPOTensor);
            }
        }    
 
        MPOOcc.orthogonalize();
        ITensor temp = TLS_Liouville::apply_MPO_for_readout(MPOOcc);
        double occ = temp.real();
        //occ /= pow(2,ttotal);
        //occ = pow(occ,0.5);
        return occ;
    }
    
    ITensor TLS_Liouville::contract_psi(MPS psi, int ttotal)
    {
        ITensor temp = psi.A(1);
        for (int j = 2; j<=(ttotal+1); j++)
        {
            temp *= psi.A(j);
        }
        return temp;
    
    }
    
    void TLS_Liouville::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::string date = TLS_Liouville::get_date_as_string();
        std::ostringstream oss1, oss2, oss3;
        oss1 << "out/timeevo_" << date.c_str() << "_Om" << Omega << "_GamTLS"<< Gamma << "_GamOut"<< GammaOut << "_GamIn" << GammaIn << "_Kappa" << Kappa << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phi" << (phi/3.141592654) << ".dat";
        file_identifier_timeevo = oss1.str();
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "w");
        fprintf(datastream_timeevo,"#time \t\t singl_ex_occ \t\t singl_gr_occ \t\t (singl_gr_occ+singl_ex_occ) \t\t trace \t\t purity \t\t contr_trace \t\t contr_ex_occ \t\t contr_gr_occ \t\t (contr_gr_occ + contr_ex_occ) \n");
        fclose(datastream_timeevo);
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
        psi = MPS(ttotal+1);
        
        // system may be in ground or exited state, timebins are initialized in vacuum
        ITensor MPSTensor = ITensor();
        for (int i=0;i<(ttotal+1);i++)
        {
            MPSTensor = ITensor(bin[i]);
            if(i==0)
            {
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
            }
            else
            {
                MPSTensor.set(bin[i](1),1.0);
            }
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
        // orthogonalize in order to establish links
        // orthocenter is now at psi.A(1)
        psi.orthogonalize();
    }
    
    void TLS_Liouville::set_up_Uevo_without_feedback()
    {
        println("set up MPO without feedback\n");
        
        Lsys = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));
        for (int p=1; p<5; p++)
        {
            Lsys.set(bin[0](1),bin[1](p),prime(bin[0](3)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](2),bin[1](p),prime(bin[0](4)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](3),bin[1](p),prime(bin[0](1)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](4),bin[1](p),prime(bin[0](2)),prime(bin[1](p)),-Cplx_i*Omega*dt);
            Lsys.set(bin[0](2),bin[1](p),prime(bin[0](1)),prime(bin[1](p)),+Cplx_i*Omega*dt);
            Lsys.set(bin[0](4),bin[1](p),prime(bin[0](3)),prime(bin[1](p)),+Cplx_i*Omega*dt);
            Lsys.set(bin[0](1),bin[1](p),prime(bin[0](2)),prime(bin[1](p)),+Cplx_i*Omega*dt);
            Lsys.set(bin[0](3),bin[1](p),prime(bin[0](4)),prime(bin[1](p)),+Cplx_i*Omega*dt);
        }
        
        Lint = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));
        
//         // Absteiger im Bad
//         Lint.set(bin[0](1),bin[1](3),prime(bin[0](3)),prime(bin[1](1)), Gamma*pow(dt,0.5));
//         Lint.set(bin[0](1),bin[1](4),prime(bin[0](3)),prime(bin[1](2)), Gamma*pow(dt,0.5));
//         Lint.set(bin[0](2),bin[1](3),prime(bin[0](4)),prime(bin[1](1)), Gamma*pow(dt,0.5));
//         Lint.set(bin[0](2),bin[1](4),prime(bin[0](4)),prime(bin[1](2)), Gamma*pow(dt,0.5));
//         
//         //Aufsteiger im Bad
//         Lint.set(bin[0](3),bin[1](1),prime(bin[0](1)),prime(bin[1](3)), -Gamma*pow(dt,0.5));
//         Lint.set(bin[0](3),bin[1](2),prime(bin[0](1)),prime(bin[1](4)), -Gamma*pow(dt,0.5));
//         Lint.set(bin[0](4),bin[1](1),prime(bin[0](2)),prime(bin[1](3)), -Gamma*pow(dt,0.5));
//         Lint.set(bin[0](4),bin[1](2),prime(bin[0](2)),prime(bin[1](4)), -Gamma*pow(dt,0.5));
//         
//         // Absteiger im Bad
//         Lint.set(bin[0](2),bin[1](1),prime(bin[0](1)),prime(bin[1](2)), -Gamma*pow(dt,0.5));
//         Lint.set(bin[0](2),bin[1](3),prime(bin[0](1)),prime(bin[1](4)), -Gamma*pow(dt,0.5));
//         Lint.set(bin[0](4),bin[1](1),prime(bin[0](3)),prime(bin[1](2)), -Gamma*pow(dt,0.5));
//         Lint.set(bin[0](4),bin[1](3),prime(bin[0](3)),prime(bin[1](4)), -Gamma*pow(dt,0.5));
//     
//         // Aufsteiger im Bad
//         Lint.set(bin[0](1),bin[1](2),prime(bin[0](2)),prime(bin[1](1)), Gamma*pow(dt,0.5));
//         Lint.set(bin[0](1),bin[1](4),prime(bin[0](2)),prime(bin[1](3)), Gamma*pow(dt,0.5));
//         Lint.set(bin[0](3),bin[1](2),prime(bin[0](4)),prime(bin[1](1)), Gamma*pow(dt,0.5));
//         Lint.set(bin[0](3),bin[1](2),prime(bin[0](4)),prime(bin[1](1)), Gamma*pow(dt,0.5));
//         
        
        // Absteiger im Bad
        Lint.set(bin[0](1),bin[1](3),prime(bin[0](3)),prime(bin[1](1)), -Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](1),bin[1](4),prime(bin[0](3)),prime(bin[1](2)), -Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](2),bin[1](3),prime(bin[0](4)),prime(bin[1](1)), -Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](2),bin[1](4),prime(bin[0](4)),prime(bin[1](2)), -Cplx_i*Gamma*pow(dt,0.5));
        
        //Aufsteiger im Bad
        Lint.set(bin[0](3),bin[1](1),prime(bin[0](1)),prime(bin[1](3)), -Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](3),bin[1](2),prime(bin[0](1)),prime(bin[1](4)), -Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](4),bin[1](1),prime(bin[0](2)),prime(bin[1](3)), -Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](4),bin[1](2),prime(bin[0](2)),prime(bin[1](4)), -Cplx_i*Gamma*pow(dt,0.5));
        
        // Absteiger im Bad
        Lint.set(bin[0](2),bin[1](1),prime(bin[0](1)),prime(bin[1](2)), Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](2),bin[1](3),prime(bin[0](1)),prime(bin[1](4)), Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](4),bin[1](1),prime(bin[0](3)),prime(bin[1](2)), Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](4),bin[1](3),prime(bin[0](3)),prime(bin[1](4)), Cplx_i*Gamma*pow(dt,0.5));
    
        // Aufsteiger im Bad
        Lint.set(bin[0](1),bin[1](2),prime(bin[0](2)),prime(bin[1](1)), Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](1),bin[1](4),prime(bin[0](2)),prime(bin[1](3)), Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](3),bin[1](2),prime(bin[0](4)),prime(bin[1](1)), Cplx_i*Gamma*pow(dt,0.5));
        Lint.set(bin[0](3),bin[1](2),prime(bin[0](4)),prime(bin[1](1)), Cplx_i*Gamma*pow(dt,0.5));
        
        
        
        LindbladLoss = ITensor(bin[0],prime(bin[0]));
        LindbladLoss.set(bin[0](4),prime(bin[0](1)), 2.0*GammaOut*dt);
        LindbladLoss.set(bin[0](2),prime(bin[0](2)), -1.0*GammaOut*dt);
        LindbladLoss.set(bin[0](3),prime(bin[0](3)), -1.0*GammaOut*dt);
        LindbladLoss.set(bin[0](4),prime(bin[0](4)), -2.0*GammaOut*dt);
        
        LindbladPump = ITensor(bin[0],prime(bin[0]));
        LindbladPump.set(bin[0](1),prime(bin[0](1)), -2.0*GammaIn*dt);
        LindbladPump.set(bin[0](2),prime(bin[0](2)), -1.0*GammaIn*dt);
        LindbladPump.set(bin[0](3),prime(bin[0](3)), -1.0*GammaIn*dt);
        LindbladPump.set(bin[0](1),prime(bin[0](4)), 2.0*GammaIn*dt);
        
        PureDephasing = ITensor(bin[0],prime(bin[0]));
        PureDephasing.set(bin[0](2), prime(bin[0](2)), -1.0*Kappa*dt);
        PureDephasing.set(bin[0](3), prime(bin[0](3)), -1.0*Kappa*dt);
        
        Lindblad = ITensor(bin[0],prime(bin[0]));
        Lindblad = LindbladLoss + LindbladPump + PureDephasing;
        
        ITensor identity_system = TLS_Liouville::identity(bin[0]);
        ITensor identity_int = TLS_Liouville::identity(bin[1]);
        
        Lindblad *= identity_int;
        
//         ITensor U_1  = Lsys + Lint + Lindblad; 
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
        
       Uevo = identity_int*identity_system + Lsys + Lint + Lindblad + 0.5*mapprime(prime(Lint)*Lint,2,1);
       //Uevo = identity_int*identity_system + Lsys + Lint + Lindblad + 0.5*mapprime(prime(Lint)*Lint,2,1) + 0.5*mapprime(prime(Lindblad)*Lindblad,2,1);
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
            Lsys.set(bin[0](2),bin[Nfb+1](p),prime(bin[0](1)),prime(bin[Nfb+1](p)),+Cplx_i*Omega*dt);
            Lsys.set(bin[0](4),bin[Nfb+1](p),prime(bin[0](3)),prime(bin[Nfb+1](p)),+Cplx_i*Omega*dt);
            Lsys.set(bin[0](1),bin[Nfb+1](p),prime(bin[0](2)),prime(bin[Nfb+1](p)),+Cplx_i*Omega*dt);
            Lsys.set(bin[0](3),bin[Nfb+1](p),prime(bin[0](4)),prime(bin[Nfb+1](p)),+Cplx_i*Omega*dt);
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
            Lfb.set(bin[0](1),bin[1](p+2),prime(bin[0](3)),prime(bin[1](p)), Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](2),bin[1](p+2),prime(bin[0](4)),prime(bin[1](p)), Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
        }
        for (int p=3; p<=4; p++)
        {
            Lfb.set(bin[0](3),bin[1](p-2),prime(bin[0](1)),prime(bin[1](p)), -Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](4),bin[1](p-2),prime(bin[0](2)),prime(bin[1](p)), -Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
        }
        for (int p=2; p<=4; p+=2)
        {
            Lfb.set(bin[0](2),bin[1](p-1),prime(bin[0](1)),prime(bin[1](p)), -Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](4),bin[1](p-1),prime(bin[0](3)),prime(bin[1](p)), -Gamma*exp(-Cplx_i*phi)*pow(dt,0.5));
        }
        for (int p=1; p<=3; p+=2)
        {
            Lfb.set(bin[0](1),bin[1](p+1),prime(bin[0](2)),prime(bin[1](p)), Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
            Lfb.set(bin[0](3),bin[1](p+1),prime(bin[0](4)),prime(bin[1](p)), Gamma*exp(Cplx_i*phi)*pow(dt,0.5));
        }
        
        LindbladLoss = ITensor(bin[0],prime(bin[0]));
        LindbladLoss.set(bin[0](4),prime(bin[0](1)), 2.0*GammaOut*dt);
        LindbladLoss.set(bin[0](2),prime(bin[0](2)), -1.0*GammaOut*dt);
        LindbladLoss.set(bin[0](3),prime(bin[0](3)), -1.0*GammaOut*dt);
        LindbladLoss.set(bin[0](4),prime(bin[0](4)), -2.0*GammaOut*dt);
        
        LindbladPump = ITensor(bin[0],prime(bin[0]));
        LindbladPump.set(bin[0](1),prime(bin[0](1)), -2.0*GammaIn*dt);
        LindbladPump.set(bin[0](2),prime(bin[0](2)), -1.0*GammaIn*dt);
        LindbladPump.set(bin[0](3),prime(bin[0](3)), -1.0*GammaIn*dt);
        LindbladPump.set(bin[0](1),prime(bin[0](4)), 2.0*GammaIn*dt);

        
        PureDephasing = ITensor(bin[0],prime(bin[0]));
        PureDephasing.set(bin[0](2), prime(bin[0](2)), -1.0*Kappa*dt);
        PureDephasing.set(bin[0](3), prime(bin[0](3)), -1.0*Kappa*dt);
        
        Lindblad = ITensor(bin[0],prime(bin[0]));
        Lindblad = LindbladLoss + LindbladPump + PureDephasing;
        
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
        TLS_Liouville::swap_feedbackbin_to_system(timestep);
        
        ITensor U,S,V,temp;
        //shift ortho-center from feedbackbin to systembin to compute observable
        temp = psi.A(timestep)*psi.A(timestep+1);
        U = ITensor(bin[0],rightLinkInd(psi,timestep+1));
        svd(temp,U,S,V,args);
        psi.setA(timestep,V);
        psi.setA(timestep+1,U*S);
        fprintf(datastream_timeevo, "%f \t %e \t %e \t %e \t %e \n", (timestep-Nfb)*dt, TLS_Liouville::excited_occupation(psi.A(timestep+1), Site), TLS_Liouville::ground_occupation(psi.A(timestep+1), Site), TLS_Liouville::trace(psi.A(timestep+1), Site), TLS_Liouville::rhosquare(psi.A(timestep+1), Site));
        if ((timestep%10)==0) 
        {
            printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", (timestep-Nfb)*dt, t_end*dt, TLS_Liouville::excited_occupation(psi.A(timestep+1), Site),  TLS_Liouville::trace(psi.A(timestep+1), Site) ,TLS_Liouville::rhosquare(psi.A(timestep+1), Site));
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
    
    void TLS_Liouville::apply_Uevo_without_feedback(int timestep)
    {
        ITensor U,S,V,temp;    
        temp = noprime(Uevo*psi.A(timestep+1)*psi.A(timestep+2));
        U = ITensor(bin[timestep+1],commonIndex(psi.A(timestep+1),psi.A(timestep),Link));
        V = ITensor(bin[0],commonIndex(psi.A(timestep+2),psi.A(timestep+3)));
        svd(temp,U,S,V,args);
        psi.setA(timestep+1,U);
        psi.setA(timestep+2,V*S);
    }
    
    void TLS_Liouville::timeevolve_system_without_feedback()
    {
        TLS_Liouville::setup_dir_and_filename();
        TLS_Liouville::setup_indexbin();
        
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "a");
        TLS_Liouville::initialize_MPS_without_feedback();
        TLS_Liouville::set_up_Uevo_without_feedback();
        
        println("start timeevolution of system without feedback\n");
        
        // Versuch: contract psi
        ITensor psi_contracted = TLS_Liouville::contract_psi(psi, ttotal);
        
        fprintf(datastream_timeevo, "%f \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", 0.0*dt, TLS_Liouville::excited_occupation(psi.A(1), Site), TLS_Liouville::ground_occupation(psi.A(1), Site), (TLS_Liouville::excited_occupation(psi.A(1), Site)+TLS_Liouville::ground_occupation(psi.A(1), Site)), TLS_Liouville::trace(psi.A(1), Site), TLS_Liouville::rhosquare(psi.A(1), Site), TLS_Liouville::trace_contracted(psi_contracted), TLS_Liouville::tr_ex_occ(psi_contracted, bin[0]), TLS_Liouville::tr_gr_occ(psi_contracted, bin[0]), (TLS_Liouville::tr_ex_occ(psi_contracted, bin[0]) + TLS_Liouville::tr_gr_occ(psi_contracted, bin[0]) ));

        for(int j = 0; j<ttotal; j++)
        {
            TLS_Liouville::apply_Uevo_without_feedback(j); // nach dieser Schleife ist das OC auf psi.A(j+2). Letzter Durchgang ist bei ttotal-1, also ist es dann auf ttotal+1, also auf system.
            psi_contracted = TLS_Liouville::contract_psi(psi, ttotal);
            PrintData(psi_contracted);
            if(j<(ttotal-1))
            {
                Uevo = Uevo*delta(bin[j+1],bin[j+2])*delta(prime(bin[j+1]),prime(bin[j+2]));
            }
            fprintf(datastream_timeevo, "%f \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", (j+1)*dt, TLS_Liouville::excited_occupation(psi.A(j+2), Site), TLS_Liouville::ground_occupation(psi.A(j+2), Site), (TLS_Liouville::excited_occupation(psi.A(j+2), Site)+TLS_Liouville::ground_occupation(psi.A(j+2), Site)), TLS_Liouville::trace(psi.A(j+2), Site), TLS_Liouville::rhosquare(psi.A(j+2), Site), TLS_Liouville::trace_contracted(psi_contracted), TLS_Liouville::tr_ex_occ(psi_contracted, bin[0]), TLS_Liouville::tr_gr_occ(psi_contracted, bin[0]), (TLS_Liouville::tr_ex_occ(psi_contracted, bin[0]) + TLS_Liouville::tr_gr_occ(psi_contracted, bin[0])));
            if ((j%10)==0) 
            {
                printf("time %.3f -- ex occ=%.3f -- gr occ=%.3f -- (ex occ+gr occ=%.3f)  -- trace=%.3f -- purity=%.3f \n", (j+1)*dt, TLS_Liouville::excited_occupation(psi.A(j+2), Site), TLS_Liouville::ground_occupation(psi.A(j+2), Site), (TLS_Liouville::excited_occupation(psi.A(j+2), Site)+TLS_Liouville::ground_occupation(psi.A(j+2), Site)), TLS_Liouville::trace(psi.A(j+2), Site), TLS_Liouville::rhosquare(psi.A(j+2), Site));
            }
        }
        fclose(datastream_timeevo);
    }
    
    void TLS_Liouville::timeevolve_system_with_feedback()
    {
        println("\n############################################################## \n \n");
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "a");
        TLS_Liouville::setup_indexbin();
        TLS_Liouville::initialize_MPS_with_feedback();
        TLS_Liouville::set_up_Uevo_with_feedback();
        println("starting timeevolution of system... \n");
        printf("time %.3f of %.3f -- occ=%.10f -- trace=%.10f -- purity=%.10f \n", 0.0*dt, t_end*dt, TLS_Liouville::excited_occupation(psi.A(1), Site),  TLS_Liouville::trace(psi.A(1), Site) ,TLS_Liouville::rhosquare(psi.A(1), Site));
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

}//eof namespace




