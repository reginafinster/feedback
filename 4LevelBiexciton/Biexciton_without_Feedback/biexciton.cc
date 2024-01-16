#include "biexciton.h"


namespace itensor{
    
    Biexciton::Biexciton(std::string infile)
    {
        // params for system and interaction
        InputGroup input = InputGroup(infile,"input");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        dt = input.getReal("time_step");
        t_end = input.getReal("Ntime");
        Omega_H = input.getReal("Omega_H");
        Delta = input.getReal("Delta");
        Gamma = input.getReal("Gamma");
        Nbin = input.getInt("Nbin");
        init_GG = input.getReal("init_GG");
        init_VV = input.getReal("init_VV");
        init_HH = input.getReal("init_HH");
        init_BB = input.getReal("init_BB");
        print_every_nth_step = input.getInt("print_every_nth_step");
        save_every_nth_step = input.getInt("save_every_nth_step");
        
        // params for spectrum
        spectrum_steps = input.getInt("spectrum_steps");
        spectrum_interval = input.getReal("spectrum_interval");   
        
        // params for feedback
        phi = input.getReal("feedbackphase_in_units_of_phi");
        Nfb = input.getInt("Nfb");

        // calculate some params
        ttotal = t_end;
        ttotal += Nfb;
        args = Args("Cutoff=",cutoff,"Maxm=",maxm);
        phi *= 3.141592654;
        
        // initialize datastream
        Biexciton::setup_dir_and_filename();
    }
    
    std::string Biexciton::get_date_as_string()
    {
        time_t t = time(0); 
        struct tm *now = localtime(&t);
        std::stringstream datestring;
        datestring << now->tm_mday << '-' << (now->tm_mon + 1) << '-' << (now->tm_year + 1900);
        return datestring.str();
    }
    
    ITensor Biexciton::identity(Index s)
    {
        ITensor my_identity = ITensor(s,prime(s));
        for(int i=1; i<=s.m(); i++)
        {
            my_identity.set(s(i),prime(s)(i),1.0);
        }
        return my_identity;
    }

    double Biexciton::occupation_system(ITensor psitensor, int system_level)
    {
        
        Index s = findtype(psitensor, Site);
        if(s.rawname()!="s")
        {
            println("You tried to read the system occupation from a bath bin. Exiting.");
            throw std::exception();
        }
        ITensor SpSm = ITensor(s, prime(s));
        SpSm.set(s(system_level), prime(s(system_level)), 1.0);
        ITensor temp = dag(prime(psitensor, s)) * (SpSm*psitensor);
        double occ = temp.real();
        return occ;
    }
    
    double Biexciton::sum_probabilities(ITensor psitensor)
    {
        double sum = Biexciton::occupation_system(psitensor,1);
        sum += Biexciton::occupation_system(psitensor,2);
        sum += Biexciton::occupation_system(psitensor,3);
        sum += Biexciton::occupation_system(psitensor,4);
        return sum;
    } 
    
    ITensor Biexciton::annihilator(Index s)
    {
        ITensor Sm = ITensor(s, prime(s));
        for(int j=1;j<Nbin;j++)
        {
            Sm.set(s(j+1), prime(s(j)), pow(j,0.5));
        }
        return Sm;
    }
    
    std::complex<double> Biexciton::calc_g1_bathA(Index timebinA, Index timebinB, ITensor temp)
    {
        ITensor Sm_temp = Biexciton::annihilator(timebinA);
        Sm_temp *= Biexciton::identity(timebinB);
        ITensor temp_tensor = dag(noprime(temp*Sm_RefBinA)) * noprime(Sm_temp*temp);
        std::complex<double> g1_cplx = temp_tensor.cplx();
        return g1_cplx;
    }
    
    std::complex<double> Biexciton::calc_g1_bathB(Index timebinB, Index timebinA, ITensor temp)
    {
        ITensor Sm_temp = Biexciton::annihilator(timebinB);
        Sm_temp *= Biexciton::identity(timebinA);
        ITensor temp_tensor = dag(noprime(temp*Sm_RefBinB)) * noprime(Sm_temp*temp);
        std::complex<double> g1_cplx = temp_tensor.cplx();
        return g1_cplx;
    }
    
    std::complex<double> Biexciton::calc_g2_bathA(Index timebinA, Index timebinB, ITensor temp)
    {
        ITensor Sm_temp = Biexciton::annihilator(timebinA);
        Sm_temp *= Biexciton::identity(timebinB); 
        ITensor temp_tensor = dag(noprime(Sm_temp*noprime(Sm_RefBinA*temp))) * noprime(Sm_temp*noprime(Sm_RefBinA*temp));
        std::complex<double> g2_cplx = temp_tensor.cplx();
        return g2_cplx;
    }
    
    std::complex<double> Biexciton::calc_g2_bathB(Index timebinB, Index timebinA, ITensor temp)
    {
        ITensor Sm_temp = Biexciton::annihilator(timebinB);
        Sm_temp *= Biexciton::identity(timebinA);
        ITensor temp_tensor = dag(noprime(Sm_temp*noprime(Sm_RefBinB*temp))) * noprime(Sm_temp*noprime(Sm_RefBinB*temp));
        std::complex<double> g2_cplx = temp_tensor.cplx();
        return g2_cplx;
    }
    
    std::complex<double> Biexciton::calc_pop_bathA(ITensor temp)
    {
        ITensor temp_tensor = dag(noprime(Sm_RefBinA*temp)) * noprime(Sm_RefBinA*temp);
        std::complex<double> pop_bath_cplx = temp_tensor.cplx();
        return pop_bath_cplx;
    }
    
    std::complex<double> Biexciton::calc_pop_bathB(ITensor temp)
    {
        ITensor temp_tensor = dag(noprime(Sm_RefBinB*temp)) * noprime(Sm_RefBinB*temp);            
        std::complex<double> pop_bath_cplx = temp_tensor.cplx();
        return pop_bath_cplx;
    }

    double Biexciton::norm(ITensor psi)
    {
        return (dag(psi)*psi).real();
    }

    void Biexciton::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        std::string date = Biexciton::get_date_as_string();
        
        std::ostringstream oss1, oss2, oss3, oss4, oss5;
        
        oss1 << "out/TLS_biexciton_timeevo" << date.c_str() << "_Om" << Omega_H << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_timeevo = oss1.str();
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "w");
        fprintf(datastream_timeevo,"#time \t popG \t popV \t popH \t popB \t sum of probab \t norm \n");
        fclose(datastream_timeevo);
        
        oss2 << "out/correlation_biexciton_V_pol" << date.c_str() << "_Om" << Omega_H << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_correlationBforV = oss2.str();
        datastream_correlationBforV = fopen(file_identifier_correlationBforV.c_str(), "w");
        fprintf(datastream_correlationBforV, "#time \t Re[g1_V] \t Im[g1_V] \t Re[g2_V] \t Im[g2_V] \t pop_V_bath^2 \n");
        fclose(datastream_correlationBforV);
        
        oss3 << "out/correlation_biexciton_H_pol" << date.c_str() << "_Om" << Omega_H << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_correlationAforH = oss3.str();
        datastream_correlationAforH = fopen(file_identifier_correlationAforH.c_str(), "w");
        fprintf(datastream_correlationAforH, "#time \t Re[g1_H] \t Im[g1_H] \t Re[g2_H] \t Im[g2_H] \t pop_V_bath^2 \n");
        fclose(datastream_correlationAforH);
        
        oss4 << "out/spectrum_biexciton_V_pol" << date.c_str() << "_Om" << Omega_H << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_spectrumBforV = oss4.str();
        datastream_spectrumBforV = fopen(file_identifier_spectrumBforV.c_str(), "w");
        fprintf(datastream_spectrumBforV, "#omega \t Re_V[S] \t Im_V[S] \t bath_pop_V \n");
        fclose(datastream_spectrumBforV);
        
        oss5 << "out/spectrum_biexciton_H_pol" << date.c_str() << "_Om" << Omega_H << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_spectrumAforH = oss5.str();
        datastream_spectrumAforH = fopen(file_identifier_spectrumAforH.c_str(), "w");
        fprintf(datastream_spectrumAforH, "#omega \t Re_V[S] \t Im_V[S] \t bath_pop_V \n");
        fclose(datastream_spectrumAforH);
        
    }
    
    void Biexciton::setup_indexbins()
    {
        // binA is vor H-polarization, binB for V-polarization
        binA = std::vector<Index>(ttotal+1);
        binB = std::vector<Index>(ttotal+1);

        // first entry of both vectors is the system
        binA.at(0) = Index("s",4,Site);
        binB.at(0) = Index("s",4,Site);
        
        // rest is filled with indices for timesteps
        for (int i=1;i<(ttotal+1);i++)
        {
            binA.at(i) = Index(nameint("kA",i),Nbin,Site);
            binB.at(i) = Index(nameint("kB",i),Nbin,Site);
        }
    }
    
    void Biexciton::initialize_MPS()
    {
        // We need one bath tensor for each timestep and one for the system
        psi = MPS(ttotal+1);
        
        // timebins are initialized in vacuum
        ITensor MPSTensor = ITensor();

        for (int i=1;i<Nfb+1;i++)
        {
            MPSTensor = ITensor(binA[i],binB[i]);
            MPSTensor.set(binA[i](1), binB[i](1), 1.0);
            psi.setA((i),MPSTensor);
        }
        
        // system is at position A(Nfb+1) initially
        MPSTensor = ITensor(binA[0]);
        MPSTensor.set(binA[0](1),init_GG);
        MPSTensor.set(binA[0](2),init_VV);
        MPSTensor.set(binA[0](3),init_HH);
        MPSTensor.set(binA[0](4),init_BB);
        printf("Initialize system with %.3f in |G>, %.3f in |V>, %.3f in |H> and %.3f in |B>\n \n", init_GG, init_VV, init_HH, init_BB);
        psi.setA((Nfb+1),MPSTensor);
        
        // rest is again filled with timebins in vacuum
        for (int i=Nfb+2;i<ttotal+2;i++)
        {
            MPSTensor = ITensor(binA[i-1],binB[i-1]);
            MPSTensor.set(binA[i-1](1), binB[i-1](1), 1.0);
            psi.setA((i),MPSTensor);
        }
        psi.orthogonalize();
    }
    
    void Biexciton::set_up_MPO()
    {
        Hsys = ITensor(binA[0],prime(binA[0]));
        Hsys.set(binA[0](2), prime(binA[0](2)), Cplx_i*Delta*pow(dt,0.5));
        Hsys.set(binA[0](3), prime(binA[0](3)), Cplx_i*Delta*pow(dt,0.5));
        
        Hsys.set(binA[0](3), prime(binA[0](1)), Cplx_i*Omega_H*pow(dt,0.5));
        Hsys.set(binA[0](1), prime(binA[0](3)), Cplx_i*Omega_H*pow(dt,0.5));
        Hsys.set(binA[0](3), prime(binA[0](4)), Cplx_i*Omega_H*pow(dt,0.5));
        Hsys.set(binA[0](4), prime(binA[0](3)), Cplx_i*Omega_H*pow(dt,0.5));
        
        // WW mit H-Polarisation über Bad A:
        HintA = ITensor(binA[0],prime(binA[0]),binA[Nfb+1],prime(binA[Nfb+1]));
        
        // WW mit V-Polarisation über Bad B:
        HintB = ITensor(binA[0],prime(binA[0]),binB[Nfb+1],prime(binB[Nfb+1]));
        //HfbB = ITensor(binB[0],prime(binB[0]),binB[1],prime(binB[1]));

        for(int i=1;i<Nbin;i++)
        {            
            // Absteiger im System, Aufsteiger im Bad B, v-polarised
            HintB.set(binA[0](2),prime(binA[0](1)),binB[Nfb+1](i),prime(binB[Nfb+1](i+1)), Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5)); 
            HintB.set(binA[0](4),prime(binA[0](2)),binB[Nfb+1](i),prime(binB[Nfb+1](i+1)), Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5)); 
            
            // Aufsteiger im System, Absteiger im Bad B, v-polarised
            HintB.set(binA[0](1),prime(binA[0](2)),binB[Nfb+1](i+1),prime(binB[Nfb+1](i)), (+1.0)*Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5)); 
            HintB.set(binA[0](2),prime(binA[0](4)),binB[Nfb+1](i+1),prime(binB[Nfb+1](i)), (+1.0)*Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5)); 
            // Hfb.set(bin[0](1),prime(bin[0](2)),bin[1](i+1),prime(bin[1](i)), pow(dt,0.5)*Gamma*exp(Cplx_i*phi)*pow(i,0.5));
            
            // Absteiger im System, Aufsteiger im Bad A, h-polarised
            HintA.set(binA[0](3),prime(binA[0](1)),binA[Nfb+1](i),prime(binA[Nfb+1](i+1)), Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5)); 
            HintA.set(binA[0](4),prime(binA[0](3)),binA[Nfb+1](i),prime(binA[Nfb+1](i+1)), Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5)); 
            // Hfb.set(bin[0](2),prime(bin[0](1)),bin[1](i),prime(bin[1](i+1)), (-1.0)*pow(dt,0.5)*Gamma*exp(-Cplx_i*phi)*pow(i,0.5));
            
            // Aufsteiger im System, Absteiger im Bad A, h-polarised
            HintA.set(binA[0](1),prime(binA[0](3)),binA[Nfb+1](i+1),prime(binA[Nfb+1](i)), (+1.0)*Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5)); 
            HintA.set(binA[0](3),prime(binA[0](4)),binA[Nfb+1](i+1),prime(binA[Nfb+1](i)), (+1.0)*Cplx_i*pow(dt,0.5)*Gamma*pow(i,0.5));
           
        }
        
        ITensor identity_system = Biexciton::identity(binA[0]);
        ITensor identity_intA = Biexciton::identity(binA[Nfb+1]);
        ITensor identity_intB = Biexciton::identity(binB[Nfb+1]);
        ITensor identity_fbA = Biexciton::identity(binA[1]);
        ITensor identity_fbB = Biexciton::identity(binB[1]);
        
        Hsys *= (identity_intA*identity_intB*identity_fbA*identity_fbB);
        HintA *= (identity_intB*identity_fbA*identity_fbB);
        HintB *= (identity_intA*identity_fbA*identity_fbB);
        
        Hint = HintA + HintB;
        
        // workaround: set up feedback operator as identity
//         Hfb = Biexciton::identity(binB[1]);
//         Hfb *= Cplx_i*pow(dt,0.5);
//         Hfb *= (identity_system*identity_fbA*identity_intA*identity_intB);
//         Hint += Hfb;
        
        ITensor U_fb_1  = Hsys + Hint; 
        ITensor U_fb_2  = (1./2.)  * mapprime(U_fb_1*prime(U_fb_1),2,1);
        ITensor U_fb_3  = (1./3.)  * mapprime(U_fb_1*prime(U_fb_2),2,1);
        ITensor U_fb_4  = (1./4.)  * mapprime(U_fb_1*prime(U_fb_3),2,1);
        ITensor U_fb_5  = (1./5.)  * mapprime(U_fb_1*prime(U_fb_4),2,1);
        ITensor U_fb_6  = (1./6.)  * mapprime(U_fb_1*prime(U_fb_5),2,1);
        ITensor U_fb_7  = (1./7.)  * mapprime(U_fb_1*prime(U_fb_6),2,1);
        ITensor U_fb_8  = (1./8.)  * mapprime(U_fb_1*prime(U_fb_7),2,1);
        ITensor U_fb_9  = (1./9.)  * mapprime(U_fb_1*prime(U_fb_8),2,1);
        ITensor U_fb_10 = (1./10.) * mapprime(U_fb_1*prime(U_fb_9),2,1);
    
        Uevo = identity_system*identity_intA*identity_intB*identity_fbA*identity_fbB + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7 + U_fb_8 + U_fb_9 + U_fb_10 ;
        
        // second order evolution
        // Hint += Hfb;
        // Uevo = identity_system*identity_int*identity_fb + Hsys + Hint + 0.5*mapprime(prime(Hint)*Hint,2,1);
    }
    
    void Biexciton::timeevolve_system()     
    {
        println("\n############################################################## \n \n");
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "a");
        Biexciton::setup_indexbins();
        Biexciton::initialize_MPS();
        println("starting timeevolution of system... \n");
        printf("time %.3f of %.3f -- norm = %.10f -- occGG = %.10f \n", 0.0*dt, (ttotal-Nfb)*dt, (1.0 - Biexciton::norm(psi.A(Nfb+1))), Biexciton::occupation_system(psi.A(Nfb+1), 1));
        Biexciton::set_up_MPO();
        for(int j=Nfb; j<ttotal; j++)
        {
            Biexciton::apply_Uevo(j);
            if(j<(ttotal-1))
            {
                // Uevo needs to act on future bin of MPS for next timestep:
                Uevo *= delta(binA[j+1],binA[j+2]);
                Uevo *= delta(binB[j+1],binB[j+2]);
                Uevo *= delta(prime(binA[j+1]),prime(binA[j+2]));
                Uevo *= delta(prime(binB[j+1]),prime(binB[j+2]));
                Uevo *= delta(binA[j-Nfb+1],binA[j-Nfb+2]);
                Uevo *= delta(prime(binA[j-Nfb+1]),prime(binA[j-Nfb+2]));
                Uevo *= delta(binB[j-Nfb+1],binB[j-Nfb+2]);
                Uevo *= delta(prime(binB[j-Nfb+1]),prime(binB[j-Nfb+2]));
            }
        }
        // system is now at the end of MPS, thus at psi.A(ttotal+1).
        // last bin with feedback interaction has index bin[Nfb] and is at position psi.A(Nfb).
        fclose(datastream_timeevo);
        println("\nfinished timeevolution of system.\n \n");
    }
    
    void Biexciton::apply_Uevo(int timestep)
    {
        Biexciton::swap_feedbackbin_to_system(timestep);
        // now feedbackbin is left of system 

        ITensor U,S,V,temp;
        //shift ortho-center from feedbackbin to systembin to compute observable
        temp = psi.A(timestep)*psi.A(timestep+1);
        U = ITensor(binA[0],rightLinkInd(psi,timestep+1));
        svd(temp,U,S,V,args);
        psi.setA(timestep,V);
        psi.setA(timestep+1,U*S);
        // system is now at A(timestep+1)
        
        if((timestep%save_every_nth_step)==0) 
        {
            fprintf(datastream_timeevo, "%f \t %e \t %e \t %e \t %e \t %e \t %e \n", (timestep-Nfb)*dt, Biexciton::occupation_system(psi.A(timestep+1), 1), Biexciton::occupation_system(psi.A(timestep+1), 2), Biexciton::occupation_system(psi.A(timestep+1), 3), Biexciton::occupation_system(psi.A(timestep+1), 4), Biexciton::sum_probabilities(psi.A(timestep+1)), (1.0 - Biexciton::norm(psi.A(timestep+1))));
        }
        
        if((timestep%print_every_nth_step)==0) 
        {
            printf("time %.3f of %.3f -- norm = %e -- occGG = %e \n",(timestep-Nfb)*dt, (ttotal-Nfb)*dt, (1.0 - Biexciton::norm(psi.A(timestep+1))), Biexciton::occupation_system(psi.A(timestep+1), 1));
        }
            
        //Index currenttimebin = findtype(psi.A(timestep+2),Site);
        // in this step, system is brought to A(timestep+2), one to the right
        // here, we look for indices of next future timebin. Because system is in MPS on the left of the timebin now, indices of bin are one less than number of tensor in MPS!
        IndexSet inds = psi.A(timestep+2).inds();
        if(findindex(inds, binA[timestep+1])<0)
        {
            println("Mismatch of indices during timeevolution, while bringing system to the right. Exiting.");
            throw std::exception();
        }
        Index currenttimebinA = inds[findindex(inds, binA[timestep+1])];
        Index currenttimebinB = inds[findindex(inds, binB[timestep+1])];
        
        temp = noprime(psi.A(timestep)*psi.A(timestep+1)*psi.A(timestep+2)*Uevo);
        
        // temp has three site-indices, thus we need to decompose twice
        // bring system to A(timestep+2), thus to the right
        U = ITensor(binA[0],leftLinkInd(psi,timestep+3));
        svd(temp,U,S,V,args);
        V *= S;
        psi.setA(timestep+2,U);
            
        // bring current timebin one position to the left (which is middle position)
        // and feedbacktimebin back to its position at A(timestep) 
        U = ITensor(currenttimebinA, currenttimebinB, commonIndex(U,V));
        temp=V;
        svd(temp,U,S,V,args);
        psi.setA(timestep+1,U);
        psi.setA(timestep,V*S);

        Biexciton::swap_feedbackbin_to_origin(timestep);
    }
    
    void Biexciton::swap_feedbackbin_to_system(int timestep)
    {
        ITensor U,S,V,temp;
        IndexSet inds = psi.A(timestep+1-Nfb).inds();
        if(findindex(inds, binA[timestep+1-Nfb])<0)
        {
            println("Mismatch of indices during timeevolution, while swapping feedbackbin to system. Exiting.");
            throw std::exception();
        }
        Index feedbacktimebinA =  inds[findindex(inds, binA[timestep+1-Nfb])];
        Index feedbacktimebinB =  inds[findindex(inds, binB[timestep+1-Nfb])];
        // Index feedbacktimebin = findtype(psi.A(timestep+1-Nfb),Site);
        for(int i=timestep+1-Nfb; i<timestep; i++)
        {
            temp = psi.A(i)*psi.A(i+1);
            U = ITensor(feedbacktimebinA, feedbacktimebinB, commonIndex(psi.A(i+1),psi.A(i+2)));
            svd(temp,U,S,V,args);
            psi.setA(i,V);
            psi.setA(i+1,U*S);
        }
    }
    
    void Biexciton::swap_feedbackbin_to_origin(int timestep)
    {
        ITensor temp,U,S,V;
        // feedbackbin carries the index with a difference of Nfb to system, thus its index is (timestep+1)-Nfb (the +1 because we start counting at Nfb)
        IndexSet inds = psi.A(timestep).inds();        
        if(findindex(inds, binA[timestep-Nfb+1])<0)
        {
            println("Mismatch of indices during timeevolution, while swapping feedbackbin to origin. Exiting");
            throw std::exception();
        }
        Index feedbacktimebinA = inds[findindex(inds, binA[timestep-Nfb+1])];
        Index feedbacktimebinB = inds[findindex(inds, binB[timestep-Nfb+1])];
        // Index feedbacktimebin = findtype(psi.A(timestep),Site);
        // loop ends when feedback bin is located one bin right of its original position
        // because in the last step, the orthocenter needs to be placed
        // at the future feedback bin for the next time step
        for(int i=timestep; i>timestep-Nfb+2; i--) 
        {
            temp = psi.A(i-1)*psi.A(i);
            U = ITensor(feedbacktimebinA, feedbacktimebinB, leftLinkInd(psi,i-1));
            svd(temp,U,S,V,args);
            psi.setA(i,V);
            psi.setA(i-1,U*S);
        }
        temp = psi.A(timestep-Nfb+1)*psi.A(timestep-Nfb+2); 
        U = ITensor(feedbacktimebinA, feedbacktimebinB, leftLinkInd(psi,timestep-Nfb+1));
        svd(temp,U,S,V,args);
        psi.setA(timestep-Nfb+2,V*S);
        psi.setA(timestep-Nfb+1,U);

    }
    
    void Biexciton::calc_correlation()
    {
        println("start calculating g(1) and g(2) functions\n");
        datastream_correlationAforH = fopen(file_identifier_correlationAforH.c_str(), "a");
        datastream_correlationBforV = fopen(file_identifier_correlationBforV.c_str(), "a");
        // set up reference annihilator:
        // first find reference bins: these are the last bins
        // which have been in feedback-interaction with the system. 
        // after time loop, these bins are located at psi.A(ttotal-Nfb)
        // and the orthocenter is at psi.A(ttotal-Nfb+1)
        // so move orthocenter to reference bins first, but do not move the bins!
        // Attention, here the index of the timebin is again the same as the number of tensor in MPS because system is at right of timebin.
        ITensor U, S, V, temp;
        
        // if spectrum of feedback bath is required: referencetimebin = ttotal-Nfb.
        // if spectrum of interaction without feedback bath: referencetimebin = ttotal.
        // system is at the end of MPS after time evolution, thus at psi.A(ttotal+1)
        // thus last bath bin is at psi.A(ttotal).
        referencetimebin = ttotal-Nfb;
        
        temp = psi.A(referencetimebin)*psi.A(referencetimebin+1);
        U = ITensor(binA[referencetimebin], binB[referencetimebin], leftLinkInd(psi,referencetimebin));
        svd(temp,U,S,V,args);
        psi.setA(referencetimebin,U*S);
        psi.setA(referencetimebin+1,V);
        
        // next, set up two reference annihilators:
        Sm_RefBinA = Biexciton::annihilator(binA[referencetimebin]);
        Sm_RefBinA *= Biexciton::identity(binB[referencetimebin]);
        
        Sm_RefBinB = Biexciton::annihilator(binB[referencetimebin]);
        Sm_RefBinB *= Biexciton::identity(binA[referencetimebin]);
        
        // allocate vectors
        g1_bathA = std::vector<Complex>(spectrum_steps+2);
        g2_bathA = std::vector<Complex>(spectrum_steps+2);
        pop_bathA = std::vector<Complex>(spectrum_steps+2);
        
        g1_bathB = std::vector<Complex>(spectrum_steps+2);
        g2_bathB = std::vector<Complex>(spectrum_steps+2);
        pop_bathB = std::vector<Complex>(spectrum_steps+2);
        
        // calculate self correlation
        // TODO wie kann ich erreichen, den Code ohne Feedback so zu schreiben, dass ich keine Veränderungen bezüglich Gamma habe?
        ITensor g1_tensorA = dag(noprime(psi.A(referencetimebin)*Sm_RefBinA)) * noprime(Sm_RefBinA*psi.A(referencetimebin));
        g1_bathA[0] = g1_tensorA.cplx();
        ITensor g1_tensorB = dag(noprime(psi.A(referencetimebin)*Sm_RefBinB)) * noprime(Sm_RefBinB*psi.A(referencetimebin));
        g1_bathB[0] = g1_tensorB.cplx();
        
        ITensor g2_tensorA = dag(noprime(Sm_RefBinA*noprime(Sm_RefBinA*psi.A(referencetimebin)))) *
        noprime(Sm_RefBinA*noprime(Sm_RefBinA*psi.A(referencetimebin)));
        g2_bathA[0] = g2_tensorA.cplx();
        ITensor g2_tensorB = dag(noprime(Sm_RefBinB*noprime(Sm_RefBinB*psi.A(referencetimebin)))) *
        noprime(Sm_RefBinB*noprime(Sm_RefBinB*psi.A(referencetimebin)));
        g2_bathB[0] = g2_tensorB.cplx();
        
        ITensor pop_bath_tensorA = dag(noprime(Sm_RefBinA*psi.A(referencetimebin))) * noprime(Sm_RefBinA*psi.A(referencetimebin));
        pop_bathA[0] = pop_bath_tensorA.cplx();
        ITensor pop_bath_tensorB = dag(noprime(Sm_RefBinB*psi.A(referencetimebin))) * noprime(Sm_RefBinB*psi.A(referencetimebin));
        pop_bathB[0] = pop_bath_tensorB.cplx();
    

        for(int i=referencetimebin;i>(referencetimebin-spectrum_steps);i--)
        {
            // calculate g1, g2 and bath population
            int k = referencetimebin-i+1;
            
            // Index s = findtype(psi.A(i-1), Site);
            IndexSet inds = psi.A(i-1).inds();
            if(findindex(inds, binA[i-1])<0)
            {
                println("Mismatch of indices during calculation of correlations. Exiting.");
                throw std::exception();
            }
            Index sA = inds[findindex(inds, binA[i-1])];
            Index sB = inds[findindex(inds, binB[i-1])];
            
            temp = psi.A(i)*psi.A(i-1);
            g1_bathA[k] = Biexciton::calc_g1_bathA(sA, sB, temp);
            g2_bathA[k] = Biexciton::calc_g2_bathA(sA, sB, temp);
            pop_bathA[k] = Biexciton::calc_pop_bathA(temp);
            
            g1_bathB[k] = Biexciton::calc_g1_bathB(sB, sA, temp);
            g2_bathB[k] = Biexciton::calc_g2_bathB(sB, sA, temp);
            pop_bathB[k] = Biexciton::calc_pop_bathB(temp);
            
            if((i%save_every_nth_step)==0) 
            {
                fprintf(datastream_correlationBforV,"%e \t %e \t %e \t %e \t %e \t %e \n",dt*k,g1_bathB[k].real(),g1_bathB[k].imag(),g2_bathB[k].real(), g2_bathB[k].imag(), pop_bathB[k].real()*pop_bathB[k].real());
                fprintf(datastream_correlationAforH,"%e \t %e \t %e \t %e \t %e \t %e \n",dt*k,g1_bathA[k].real(),g1_bathA[k].imag(),g2_bathA[k].real(), g2_bathA[k].imag(), pop_bathA[k].real()*pop_bathA[k].real());
            }
            
            if((i%print_every_nth_step)==0) 
            {
                printf("g1_bath[%i]=%e+i%e -- g2_bath[%i]=%e+i%e -- (bath_pop)^2=%e \n", k, g1_bathB[k].real(), g1_bathB[k].imag(), k, g2_bathB[k].real(), g2_bathB[k].imag(), pop_bathB[k].real()*pop_bathB[k].real());
            }
            
            // swap bin[referencetimebin] and orthocenter one step to the left for next step of calculation
            if(i>(referencetimebin-spectrum_steps+1))
            {
                ITensor U1, V1, S1;
                U1 = ITensor(binA[referencetimebin], binB[referencetimebin], commonIndex(psi.A(i-1), psi.A(i-2)));
                svd(temp,U1,S1,V1);
                psi.setA(i-1,U1*S1);
                psi.setA(i,V1);
            }
      }
        fclose(datastream_correlationAforH);
        fclose(datastream_correlationBforV);
        println("\nfinished calculating g(1) and g(2) functions\n");
    }
    
    
    void Biexciton::calc_spectrum()
    {
        println("start calculating fourier trafo for spectrum....");
        // params for normalization and calculation
        spectrumA = std::vector<Complex>(spectrum_steps+2);
        spectrumB = std::vector<Complex>(spectrum_steps+2);
        dw = spectrum_interval/spectrum_steps;
        
        // needed to calculate the spectrum without coherent part
        std::complex<double> b_offsetA = g1_bathA[spectrum_steps];
        std::complex<double> b_offsetB = g1_bathB[spectrum_steps];
        // needed for normalization TODO correct would be to use (bath_pop)^2, so check
        double g2_normA = g2_bathA[spectrum_steps].real(); 
        double g2_normB = g2_bathB[spectrum_steps].real(); 
        datastream_spectrumAforH = fopen(file_identifier_spectrumAforH.c_str(), "a");
        datastream_spectrumBforV = fopen(file_identifier_spectrumBforV.c_str(), "a");

        // calculate spectral density with Wiener-Chinschin
        // thus make a fourier trafo of g(1)
        for(int i=0;i<=spectrum_steps;i++)
        {
            om = -0.5*(spectrum_interval)+i*dw;
            spectrumA[i] = 0.0;
            spectrumB[i] = 0.0;
            for(int k=0;k<=spectrum_steps;k++)
            {
                spectrumA[i] += exp(-Cplx_i*om*k*dt)*(g1_bathA[k]-b_offsetA);
                spectrumB[i] += exp(-Cplx_i*om*k*dt)*(g1_bathB[k]-b_offsetB); 

            }
            fprintf(datastream_spectrumAforH,"%e \t %e \t %e \t %e \n", om, spectrumA[i].real(), spectrumA[i].imag(), g2_bathA[i].real()/g2_normA);
            fprintf(datastream_spectrumBforV,"%e \t %e \t %e \t %e \n", om, spectrumB[i].real(), spectrumB[i].imag(), g2_bathB[i].real()/g2_normB);
        }
        fclose(datastream_spectrumAforH);
        fclose(datastream_spectrumBforV);
        println("...done\n");
    }

}//eof namespace

