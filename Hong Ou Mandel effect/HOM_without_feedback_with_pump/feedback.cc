#include "feedback.h"


namespace itensor{
    
    Feedback::Feedback(std::string infile)
    {
        InputGroup input = InputGroup(infile,"input");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        dt = input.getReal("time_step");
        t_end = input.getReal("time_end");
        Omega1 = input.getReal("Omega1");
        Gamma1 = input.getReal("Gamma1");
        Gamma2 = input.getReal("Gamma2");
        Gamma2Max = input.getReal("Gamma2Max");
        GammaSteps = input.getInt("GammaSteps");
        phi1 = input.getReal("feedbackphase1",0);
        phi1 = input.getReal("feedbackphase2",0);
        Nbin = input.getInt("Dim_of_bath");
        Nfb1 = input.getInt("feedback_time1");
        Nfb2 = input.getInt("feedback_time2");
        save_every_nth_step = input.getInt("save_every_nth_step");
        ttotal = t_end/dt;
        //ttotal += Nfb;
        args = Args("Cutoff=",cutoff,"Maxm=",maxm);
        
        phi1 *= 3.141592654;
        phi2 *= 3.141592654;

    }
    
    std::string Feedback::get_date_as_string()
    {
        time_t t = time(0); 
        struct tm *now = localtime(&t);
        std::stringstream datestring;
        datestring << now->tm_mday << '-' << (now->tm_mon + 1) << '-' << (now->tm_year + 1900);
        return datestring.str();
    }
    
    ITensor Feedback::identity(Index s)
    {
        ITensor my_identity = ITensor(s,prime(s));
        for(int i=1; i<=s.m(); i++)
        {
            my_identity.set(s(i),prime(s)(i),1.0);
        }
        return my_identity;
    }

    double Feedback::TLS_occupation(ITensor psitensor, IndexType type)
    {
        Index s = findtype(psitensor, type);
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(2), prime(s(2)), 1.0);
        ITensor temp = dag(prime(psitensor, Site)) * (SpSm*psitensor);
        double occ = temp.real();
        return occ;
    }

    double Feedback::calc_bath_pop(ITensor temp, ITensor Sm_RefBin)
    {
        ITensor pop_bath_tensor = dag(noprime(Sm_RefBin*temp)) * noprime(Sm_RefBin*temp);            
        double pop_bath = pop_bath_tensor.real();
        return pop_bath;
    }
    
    double Feedback::norm(ITensor psi, IndexType type)
    {
        return (dag(psi)*psi).real();
    }

    ITensor Feedback::annihilator(Index s)
    {
        ITensor Sm = ITensor(s, prime(s));
        for(int j=1;j<Nbin;j++)
        {
            Sm.set(s(j+1), prime(s(j)), pow(j,0.5));
        }
        return Sm;
    }
    
    ITensor Feedback::creator(Index s)
    {
        ITensor Sm = ITensor(s, prime(s));
        for(int j=1;j<Nbin;j++)
        {
            Sm.set(s(j), prime(s(j+1)), pow(j,0.5));
        }
        return Sm;
    }
    
    void Feedback::setup_dir_and_filename()
    {
        mkdir("./out1", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        std::string date = Feedback::get_date_as_string();
        
        std::ostringstream oss1;
        oss1 << "out1/TLS1_Feedback_" << date.c_str() << "_Gam"<< Gamma1 << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb1*dt) << "_phase" << (phi1/3.141592654) << ".dat";
        file_identifier1 = oss1.str();
        datastream1 = fopen(file_identifier1.c_str(), "w");
        fprintf(datastream1,"#time \tpop_system \tpop_bath \t 1-norm \n");
        fclose(datastream1);
        
        std::ostringstream oss2;
        oss2 << "out1/TLS2_Feedback_" << date.c_str() << "_Gam"<< Gamma2 << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb2*dt) << "_phase" << (phi2/3.141592654) << ".dat";
        file_identifier2 = oss2.str();
        datastream2 = fopen(file_identifier2.c_str(), "w");
        fprintf(datastream2,"#time \tpop_system \tpop_bath \t 1-norm \n");
        fclose(datastream2);
        
        std::ostringstream oss3;
        oss3 << "out1/correlation_" << date.c_str() << "_Gam1"<< Gamma1 << "_Gam2"<< Gamma2 << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime1" << (Nfb1*dt) << "_phase1" << (phi1/3.141592654) << ".dat";
        file_identifier3 = oss3.str();
        datastream3 = fopen(file_identifier3.c_str(), "w");
        fprintf(datastream3,"#tau \t g(2) of various values for td \n");
        fclose(datastream3);
        
//         std::ostringstream oss4;
//         oss4 << "out/bath_pop_" << date.c_str() << "_Gam1"<< Gamma1 << "_Gam2"<< Gamma2 << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime1" << (Nfb1*dt) << "_phase1" << (phi1/3.141592654) << ".dat";
//         file_identifier4 = oss4.str();
//         datastream4 = fopen(file_identifier4.c_str(), "w");
//         fprintf(datastream4,"#td \t pop1 \t pop2 \n");
//         fclose(datastream4);
        
//         std::ostringstream oss6;
//         oss6 << "out/g2_tau_td_fix_" << date.c_str() << "_Gam1"<< Gamma1 << "_Gam2"<< Gamma2 << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime1" << (Nfb1*dt) << "_phase1" << (phi1/3.141592654) << ".dat";
//         file_identifier6 = oss6.str();
//         datastream6 = fopen(file_identifier6.c_str(), "w");
//         fprintf(datastream6,"#tau \t sum \t t^2-part \t tr-part \t rt-part \t r^2-part \n");
//         fclose(datastream6);
        
    }
    
    void Feedback::setup_indexbins()
    {
        bin1 = std::vector<Index>(ttotal+1);
        bin2 = std::vector<Index>(ttotal+1);
        
        // first entry of vector is the system
        bin1.at(0) = Index("s1",2,Site);
        bin2.at(0) = Index("s2",2,Site);
        
        // rest is filled with indices for timesteps
        for (int i=1;i<(ttotal+1);i++)
        {
            bin1.at(i) = Index(nameint("k1_",i),Nbin,Site);
            bin2.at(i) = Index(nameint("k2_",i),Nbin,Site);
        }
    }
    
    void Feedback::initialize_MPS1_without_feedback()
    {
        // We need one tensor for each timestep and first one is system
        psi1 = MPS(ttotal+1);
        
        // system may be in ground or exited state, timebins are initialized in vacuum
        ITensor MPSTensor = ITensor();
        for (int i=0;i<(ttotal+1);i++)
        {
            MPSTensor = ITensor(bin1[i]);
            if(i==0)
            {
                MPSTensor.set(bin1[i](2),1.0);
                println("Initialize system1 in excited state\n");
            }
            else
            {
                MPSTensor.set(bin1[i](1),1.0);
            }
            psi1.setA((i+1),MPSTensor);
        }
        psi1.orthogonalize();
    }
    
    void Feedback::initialize_MPS2_without_feedback()
    {
        // We need one tensor for each timestep and first one is system
        psi2 = MPS(ttotal+1);
        
        // system may be in ground or exited state, timebins are initialized in vacuum
        ITensor MPSTensor = ITensor();
        for (int i=0;i<(ttotal+1);i++)
        {
            MPSTensor = ITensor(bin2[i]);
            if(i==0)
            {
                MPSTensor.set(bin2[i](2),1.0);
                println("Initialize system2 in excited state\n");
            }
            else
            {
                MPSTensor.set(bin2[i](1),1.0);
            }
            psi2.setA((i+1),MPSTensor);
        }
        psi2.orthogonalize();
    }
    
    void Feedback::set_up_MPO1_without_feedback()
    {
        println("set up MPO1 without feedback\n");
        
        Hsys1 = ITensor(bin1[0],prime(bin1[0]));
        Hsys1.set(bin1[0](1), prime(bin1[0](2)),-Cplx_i*Omega1*dt);
        Hsys1.set(bin1[0](2), prime(bin1[0](1)),-Cplx_i*Omega1*dt);
        
        Hint1 = ITensor(bin1[0],prime(bin1[0]),bin1[1],prime(bin1[1]));
        for (int i=1; i<Nbin; i++)
        {
            Hint1.set(bin1[0](2),prime(bin1[0](1)),bin1[1](i),prime(bin1[1](i+1)), Gamma1*pow(dt,0.5)*pow(i,0.5));
            Hint1.set(bin1[0](1),prime(bin1[0](2)),bin1[1](i+1),prime(bin1[1](i)), (-1.0)*Gamma1*pow(dt,0.5)*pow(i,0.5));
        }
        
        ITensor identity_system1 = Feedback::identity(bin1[0]);
        ITensor identity_int1 = Feedback::identity(bin1[1]);
        
        ITensor U_1_1  = Hint1 + Hsys1*identity_int1;
        ITensor U_2_1  = (1./2.)  * mapprime(U_1_1*prime(U_1_1),2,1);
        ITensor U_3_1  = (1./3.)  * mapprime(U_1_1*prime(U_2_1),2,1);
        ITensor U_4_1  = (1./4.)  * mapprime(U_1_1*prime(U_3_1),2,1);
        ITensor U_5_1  = (1./5.)  * mapprime(U_1_1*prime(U_4_1),2,1);
        ITensor U_6_1  = (1./6.)  * mapprime(U_1_1*prime(U_5_1),2,1);
        ITensor U_7_1  = (1./7.)  * mapprime(U_1_1*prime(U_6_1),2,1);
        ITensor U_8_1  = (1./8.)  * mapprime(U_1_1*prime(U_7_1),2,1);
        ITensor U_9_1  = (1./9.)  * mapprime(U_1_1*prime(U_8_1),2,1);
        ITensor U_10_1 = (1./10.) * mapprime(U_1_1*prime(U_9_1),2,1);
        
        Uevo1 = identity_int1*identity_system1 + U_1_1+U_2_1+U_3_1+U_4_1+U_5_1+U_6_1+U_7_1+U_8_1+U_9_1+U_10_1;
        //Uevo = identity_int*identity_system + Hsys*identity_int + Hint + 0.5*mapprime(prime(Hint)*Hint,2,1);
    }
    
    void Feedback::set_up_MPO2_without_feedback()
    {        
        println("set up MPO2 without feedback\n");
        Hint2 = ITensor(bin2[0],prime(bin2[0]),bin2[1],prime(bin2[1]));
        for (int i=1; i<Nbin; i++)
        {
            Hint2.set(bin2[0](2),prime(bin2[0](1)),bin2[1](i),prime(bin2[1](i+1)), Gamma2*pow(dt,0.5)*pow(i,0.5));
            Hint2.set(bin2[0](1),prime(bin2[0](2)),bin2[1](i+1),prime(bin2[1](i)), (-1.0)*Gamma2*pow(dt,0.5)*pow(i,0.5));
        }
        
        ///////////////////////// ohne FB aus Mollow-Code ////////////////////////////
        // Hint2 = ITensor(bin2[0],prime(bin2[0]),bin2[1],prime(bin2[1]));
        // for (int i=1; i<Nbin; i++)
        // {
            // Hint2.set(bin2[0](2),prime(bin2[0](1)),bin2[1](i),prime(bin2[1](i+1)), -Cplx_i*Gamma2*pow(dt,0.5)*pow(i,0.5));
            // Hint2.set(bin2[0](1),prime(bin2[0](2)),bin2[1](i+1),prime(bin2[1](i)), -Cplx_i*Gamma2*pow(dt,0.5)*pow(i,0.5));
        // }
        //////////////////////////////////////////////////////
        
        ITensor identity_system2 = Feedback::identity(bin2[0]);
        ITensor identity_int2 = Feedback::identity(bin2[1]);
        
        ITensor U_1_2  = Hint2;
        ITensor U_2_2  = (1./2.)  * mapprime(U_1_2*prime(U_1_2),2,1);
        ITensor U_3_2  = (1./3.)  * mapprime(U_1_2*prime(U_2_2),2,1);
        ITensor U_4_2  = (1./4.)  * mapprime(U_1_2*prime(U_3_2),2,1);
        ITensor U_5_2  = (1./5.)  * mapprime(U_1_2*prime(U_4_2),2,1);
        ITensor U_6_2  = (1./6.)  * mapprime(U_1_2*prime(U_5_2),2,1);
        ITensor U_7_2  = (1./7.)  * mapprime(U_1_2*prime(U_6_2),2,1);
        ITensor U_8_2  = (1./8.)  * mapprime(U_1_2*prime(U_7_2),2,1);
        ITensor U_9_2  = (1./9.)  * mapprime(U_1_2*prime(U_8_2),2,1);
        ITensor U_10_2 = (1./10.) * mapprime(U_1_2*prime(U_9_2),2,1);
        
        Uevo2 = identity_int2*identity_system2 + U_1_2+U_2_2+U_3_2+U_4_2+U_5_2+U_6_2+U_7_2+U_8_2+U_9_2+U_10_2;
        
        //Uevo = identity_int*identity_system + Hsys*identity_int + Hint + 0.5*mapprime(prime(Hint)*Hint,2,1);
    }
    
    void Feedback::timeevolve_system_without_feedback()
    // past=j as MPS starts at 1 and bin[] at 0
    // system is always at psi(j+1), future at psi(j+2)
    // Uevo acts on system and future timebin
    // thus we contract these two tensors, apply Uevo and decompose
    // with U(futurebin, link to past), V(systembin, link to future)
    // and set orthocenter to V 
    // U becomes new past bin in MPS
    // V becomes new systembin at next timestep in MPS
    {
        Feedback::setup_dir_and_filename();
        Feedback::setup_indexbins();
        
        datastream1 = fopen(file_identifier1.c_str(), "a");
        datastream2 = fopen(file_identifier2.c_str(), "a");
        
        Feedback::initialize_MPS1_without_feedback();
        Feedback::initialize_MPS2_without_feedback();
        Feedback::set_up_MPO1_without_feedback();
        Feedback::set_up_MPO2_without_feedback();
        
        println("start timeevolution of system without feedback\n");
        for(int j = 0; j<ttotal; j++)
        {
            Feedback::apply_Uevo1_without_feedback(j);
            Feedback::apply_Uevo2_without_feedback(j);
            if(j<(ttotal-1))
            {
                Uevo1 = Uevo1*delta(bin1[j+1],bin1[j+2])*delta(prime(bin1[j+1]),prime(bin1[j+2]));
                Uevo2 = Uevo2*delta(bin2[j+1],bin2[j+2])*delta(prime(bin2[j+1]),prime(bin2[j+2]));
            }
            fprintf(datastream1, "%f \t %e \t %e \t %e \n", (j+1)*dt, Feedback::TLS_occupation(psi1.A(j+2), Site), temp_bath_occupation1, (1.0 - Feedback::norm(psi1.A(j+2), Site)));
            fprintf(datastream2, "%f \t %e \t %e \t %e \n", (j+1)*dt, Feedback::TLS_occupation(psi2.A(j+2), Site), temp_bath_occupation2, (1.0 - Feedback::norm(psi2.A(j+2), Site)));

            if ((j%10)==0) 
            {
                printf("time %.3f of %.3f -- norm=%.10f \n", j*dt, t_end*1., (1.0 - Feedback::norm(psi1.A(j+2), Site)));
            }
        }
        
        fclose(datastream1);
        fclose(datastream2);
    }
    
//     void Feedback::timeevolve_system2_without_feedback()
//     {
//         datastream2 = fopen(file_identifier2.c_str(), "a");
//         Feedback::initialize_MPS2_without_feedback();
//         Feedback::set_up_MPO2_without_feedback();
//         
//         println("start timeevolution of system 2 without feedback\n");
//         for(int j = 0; j<ttotal; j++)
//         {
//             Feedback::apply_Uevo2_without_feedback(j);
//             if(j<(ttotal-1))
//             {
//                 Uevo2 = Uevo2*delta(bin2[j+1],bin2[j+2])*delta(prime(bin2[j+1]),prime(bin2[j+2]));
//             }
//             fprintf(datastream2, "%f \t %e \t %e \n", j*dt, Feedback::TLS_occupation(psi2.A(j+2), Site), (1.0 - Feedback::norm(psi2.A(j+2), Site)));
// 
//             if ((j%10)==0) 
//             {
//                 printf("system2: time %.3f of %.3f -- norm=%.10f \n", j*dt, t_end*1., (1.0 - Feedback::norm(psi1.A(j+2), Site)));
//             }
//         }
//         fclose(datastream2);
//     }
    
    void Feedback::apply_Uevo1_without_feedback(int timestep)
    {
        ITensor U,S,V,temp;    
        temp = noprime(Uevo1*psi1.A(timestep+1)*psi1.A(timestep+2));
        U = ITensor(bin1[timestep+1],commonIndex(psi1.A(timestep+1),psi1.A(timestep),Link));
        V = ITensor(bin1[0],commonIndex(psi1.A(timestep+2),psi1.A(timestep+3)));
        svd(temp,U,S,V);
        temp_bath_occupation1 = Feedback::TLS_occupation(U*S, Site);
        psi1.setA(timestep+1,U);
        psi1.setA(timestep+2,V*S);
    }
    
    void Feedback::apply_Uevo2_without_feedback(int timestep)
    {
        ITensor U,S,V,temp;    
        temp = noprime(Uevo2*psi2.A(timestep+1)*psi2.A(timestep+2));
        U = ITensor(bin2[timestep+1],commonIndex(psi2.A(timestep+1),psi2.A(timestep),Link));
        V = ITensor(bin2[0],commonIndex(psi2.A(timestep+2),psi2.A(timestep+3)));
        svd(temp,U,S,V);
        temp_bath_occupation2 = Feedback::TLS_occupation(U*S, Site);
        psi2.setA(timestep+1,U);
        psi2.setA(timestep+2,V*S);
    }
    
    void Feedback::calc_g2_without_feedback()
    {
        println("start calculating g(2) function\n");
        datastream3 = fopen(file_identifier3.c_str(), "a");
        // datastream4 = fopen(file_identifier4.c_str(), "a");
        //datastream6 = fopen(file_identifier6.c_str(), "a");
        
        // find reference bin for respective t_d: for t_d = 0, it is the first bin of the MPS, thus the bin with the first emitted signal. 
        // However, the orthocenter is at system which is at psi.A(ttotal+1).
        // so move orthocenter to beginning of the MPS first, but do not move the bins!
        
        println("shift OC and allocate vectors\n");
        int referencetimebin1 = 1;
        int orthocenter1 = ttotal+1;
        ITensor U1, S1, V1, temp1;
        for(int i=orthocenter1; i>1; i--) 
        {
            temp1 = psi1.A(i-1)*psi1.A(i);
            U1 = ITensor(findtype(psi1.A(i-1),Site), leftLinkInd(psi1,i-1));
            svd(temp1,U1,S1,V1,args);
            psi1.setA(i,V1);
            psi1.setA(i-1,U1*S1);
        }

        int referencetimebin2 = 1;
        int orthocenter2 = ttotal+1;
        ITensor U2, S2, V2, temp2;
        for(int i=orthocenter2; i>1; i--) 
        {
            temp2 = psi2.A(i-1)*psi2.A(i);
            U2 = ITensor(findtype(psi2.A(i-1),Site), leftLinkInd(psi2,i-1));
            svd(temp2,U2,S2,V2,args);
            psi2.setA(i,V2);
            psi2.setA(i-1,U2*S2);
        }
        
        // allocate vectors of length of relevant MPS-part
        // in case of no feedback, this is ttotal
        int spectrum_steps = ttotal;
        
        // matrix for all g2 values (loops over tD and tau)
        g2 = std::vector<std::vector<double>>(spectrum_steps, std::vector<double>(spectrum_steps));
        
        // vectors for reference bath pop for both systems
        pop1 = std::vector<double>(spectrum_steps);
        pop2 = std::vector<double>(spectrum_steps);
        
        printf("start loop over refbin and taubin for Gamma2 = %f\n\n", Gamma2);
        // loop over all ref_bins until end of MPS
        for(int j=referencetimebin1; j<ttotal; j++)
        {
            // calculate position in vector:
            int m = j-1;
            
            // define relationship between referencetimebin1 and referencetimebin2 for the loop
            // j is referencetimebin1, l is referencetimebin2
            int l = j;
            
            //printf("start loop over ref_bin, Index is at position %i\n\n", j);
            // set up reference annihilators for this t_d step
            Sm_RefBin1 = Feedback::annihilator(bin1[j]);
            Sm_RefBin2 = Feedback::annihilator(bin2[l]);
            Sp_RefBin1 = Feedback::creator(bin1[j]);
            Sp_RefBin2 = Feedback::creator(bin2[l]);
            
            // calculate bath pop for this reference bin
            pop1[m] = Feedback::calc_bath_pop(psi1.A(j),Sm_RefBin1);
            pop2[m] = Feedback::calc_bath_pop(psi2.A(j),Sm_RefBin2);
            
            // correlation for t_d=tau is always zero 
            g2[m][0] = 0.0;
            
            for(int i=j;i<ttotal;i++)
            {
                // printf("start loop over tau_bin, Index is at position %i\n\n", i);
                // i is position of refbin, i+1 is position of taubin
                
                // set up required operators for this tau-step:
                ITensor Sm_TauBin1 = Feedback::annihilator(bin1[i+1]);
                ITensor Sm_TauBin2 = Feedback::annihilator(bin2[i+1]);
                ITensor Sp_TauBin1 = Feedback::creator(bin1[i+1]);
                ITensor Sp_TauBin2 = Feedback::creator(bin2[i+1]);
                
                // g(2) consists of four parts:
                // 1. psi1.A(ref_bin)*psi2.A(tau_bin)
                ITensor g2_1 = noprime(Sp_RefBin1*noprime(Sp_TauBin2*noprime(Sm_TauBin2*noprime(Sm_RefBin1*psi2.A(i)*psi1.A(i+1)*psi2.A(i+1)*psi1.A(i)))));
                g2_1 *= dag(psi2.A(i)*psi1.A(i+1)*psi2.A(i+1)*psi1.A(i));
                
                // 2. psi2.A(ref_bin)*psi1.A(tau_bin)*psi2.A(tau_bin)*psi1.A(ref_bin)
                ITensor g2_2 = noprime(Sp_TauBin1*noprime(Sp_RefBin2*noprime(Sm_TauBin2*noprime(Sm_RefBin1*(psi2.A(i)*psi1.A(i+1)*psi2.A(i+1)*psi1.A(i)))))); 
                g2_2 *= dag(psi2.A(i)*psi1.A(i+1)*psi2.A(i+1)*psi1.A(i));
                
                // 3. psi1.A(ref_bin)*psi2.A(tau_bin)*psi1.A(tau_bin)*psi2.A(ref_bin)
                ITensor g2_3 = noprime(Sp_TauBin2*noprime(Sp_RefBin1*noprime(Sm_TauBin1*noprime(Sm_RefBin2*(psi1.A(i)*psi2.A(i+1)*psi1.A(i+1)*psi2.A(i)))))); 
                g2_3 *= dag(psi1.A(i)*psi2.A(i+1)*psi1.A(i+1)*psi2.A(i));
                
                // 4. psi1.A(tau_bin)*psi2.A(ref_bin)
                ITensor g2_4 = noprime(Sp_RefBin2*noprime(Sp_TauBin1*noprime(Sm_TauBin1*noprime(Sm_RefBin2*psi2.A(i)*psi1.A(i+1)*psi2.A(i+1)*psi1.A(i)))));
                g2_4 *= dag(psi2.A(i)*psi1.A(i+1)*psi2.A(i+1)*psi1.A(i));

                
                g2[m][i] = (g2_1.real() - g2_2.real() - g2_3.real() + g2_4.real());
                
//                 if(m==0)
//                 {
//                     fprintf(datastream6,"%e \t %e \t %e \t %e \t %e \t %e \n", dt*i, g2[m][i], g2_1.real(), g2_2.real(), g2_3.real(), g2_4.real());
//                 }
                
                // swap bin[referencetimebin] and orthocenter one step to the right for next step of calculation
                ITensor temp1, U1, V1, S1;
                temp1 = psi1.A(i)*psi1.A(i+1);
                U1 = ITensor(bin1[j],rightLinkInd(psi1, i+1));
                svd(temp1,U1,S1,V1);
                psi1.setA(i+1,U1*S1);
                psi1.setA(i,V1);
                
                ITensor temp2, U2, V2, S2;
                temp2 = psi2.A(i)*psi2.A(i+1);
                U2 = ITensor(bin2[l],rightLinkInd(psi2, i+1));
                svd(temp2,U2,S2,V2);
                psi2.setA(i+1,U2*S2);
                psi2.setA(i,V2);
            } // eof loop over tau_bins. ref bin and orthocenter now are at psi.A(ttotal).
            
            // now swap ref_bin back to its origin at position j
            // and in the last step, place orthocenter on the next ref_bin at j+1
            ITensor U,S,V,temp;
            
            if(j<ttotal-1)
            {
                // for MPS1
                Index refbin_index1 = findtype(psi1.A(ttotal),Site);
                for(int i=ttotal; i>(j+1); i--) 
                {
                    temp = psi1.A(i-1)*psi1.A(i);
                    U = ITensor(refbin_index1, leftLinkInd(psi1,i-1));
                    svd(temp,U,S,V,args);
                    psi1.setA(i,V);
                    psi1.setA(i-1,U*S);
                }
                temp = psi1.A(j)*psi1.A(j+1);
                U = ITensor(refbin_index1, leftLinkInd(psi1,j));
                svd(temp,U,S,V,args);
                psi1.setA(j+1,V*S);
                psi1.setA(j,U);

                // and for MPS2
                Index refbin_index2 = findtype(psi2.A(ttotal),Site);
                for(int i=ttotal; i>(j+1); i--) 
                {
                    temp = psi2.A(i-1)*psi2.A(i);
                    U = ITensor(refbin_index2, leftLinkInd(psi2,i-1));
                    svd(temp,U,S,V,args);
                    psi2.setA(i,V);
                    psi2.setA(i-1,U*S);
                }
                temp = psi2.A(j)*psi2.A(j+1);
                U = ITensor(refbin_index2, leftLinkInd(psi2,j));
                svd(temp,U,S,V,args);
                psi2.setA(j+1,V*S);
                psi2.setA(j,U);
            }
        
        } // eof loop over ref_bins
        
        println("write results to file\n");
//         // write bath pop for every td to file
//         for (int m=0; m<spectrum_steps; m++)
//         {
//             fprintf(datastream4,"%e \t", dt*m);
//             fprintf(datastream4,"%e \t %e \n", pop1[m], pop2[m]);
//         }

        // write G(2) for various td and tau to file
        for (int k=0; k<spectrum_steps; k++)
        {
            for (int m=0; m<spectrum_steps; m++)
            {
                if(m==0){fprintf(datastream3,"%e \t", dt*k);}
                if(m%save_every_nth_step==0){fprintf(datastream3,"%e \t", g2[m][k]);}
            }
            fprintf(datastream3,"\n");
        }
        
        println("and integrate over all values\n");
        // integrate over t_d and tau
        std::vector<double> g2_tau_sum = std::vector<double>(spectrum_steps);
        g2_ref_sum = 0.0;
        norm1 = 0.0;
        norm2 = 0.0;
        g2_ref_sum_norm = 0.0;
        //loop for ref_bin_sum
        for (int m=0; m<spectrum_steps; m++)
        {
            g2_tau_sum[m] = 0.0;
            // loop for tau_bin_sum
            for (int k=0; k<spectrum_steps-m; k++)
            {
                g2_tau_sum[m] += g2[m][k];
            }
            g2_ref_sum += g2_tau_sum[m];
            // calculate norm:
            norm1 += pop1[m];
            norm2 += pop2[m];
        }
        g2_ref_sum_norm = g2_ref_sum/(norm1*norm2);
        
        //fclose(datastream4);
        fclose(datastream3);
        //fclose(datastream6);
        println("\nfinished calculating g(2) function\n");
    }

    
    void Feedback::loop_over_Gamma()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        std::string date = Feedback::get_date_as_string();
        
        std::ostringstream oss5;
        oss5 << "out/g2_gamma_" << date.c_str() << "_Gam1"<< Gamma1 << "_Gam2"<< Gamma2 << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime1" << (Nfb1*dt) << "_phase1" << (phi1/3.141592654) << ".dat";
        file_identifier5 = oss5.str();
        datastream5 = fopen(file_identifier5.c_str(), "w");
        fprintf(datastream5,"#Gamma2 \t g(2) \n");
        fclose(datastream5);
        
        datastream5 = fopen(file_identifier5.c_str(), "a");
        
        double delta_gamma = (Gamma2Max-Gamma2)/GammaSteps;
        for (int i=0; i<GammaSteps+1; i++)
        {
            Feedback::timeevolve_system_without_feedback();
            Feedback::calc_g2_without_feedback();
            fprintf(datastream5,"%e \t %e \t %e \t %e \t %e \n", Gamma2, g2_ref_sum_norm, g2_ref_sum, norm1, norm2);
            Gamma2 += delta_gamma;
        }
        fclose(datastream5);
    }

}//eof namespace




