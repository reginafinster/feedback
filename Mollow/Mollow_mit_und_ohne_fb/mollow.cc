#include "mollow.h"


namespace itensor{
    
    Mollow::Mollow(std::string infile)
    {
        // params for system and interaction
        InputGroup input = InputGroup(infile,"input");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        dt = input.getReal("time_step");
        t_end = input.getReal("Ntime");
        Omega = input.getReal("Omega_TLS");
        Gamma = input.getReal("Gamma");
        Nbin = input.getInt("Nbin");
        init_22 = input.getReal("init_22");
        init_11 = input.getReal("init_11");
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
        Mollow::setup_dir_and_filename();
//         Mollow::timeevolve_system();
//         Mollow::calc_correlation();
//         Mollow::calc_spectrum();
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
    
    ITensor Mollow::annihilator(Index s)
    {
        ITensor Sm = ITensor(s, prime(s));
        for(int j=1;j<Nbin;j++)
        {
            Sm.set(s(j+1), prime(s(j)), pow(j,0.5));
        }
        return Sm;
    }
    
    std::complex<double> Mollow::calc_g1_bath(Index timebin, ITensor temp)
    {
        ITensor Sm_temp = Mollow::annihilator(timebin);
        g1_tensor = dag(noprime(temp*Sm_RefBin)) * noprime(Sm_temp*temp);
        std::complex<double> g1_cplx = g1_tensor.cplx();
        return g1_cplx;
    }
    
    std::complex<double> Mollow::calc_g2_bath(Index timebin, ITensor temp)
    {
        ITensor Sm_temp = Mollow::annihilator(timebin);
        g2_tensor = dag(noprime(Sm_temp*noprime(Sm_RefBin*temp))) * noprime(Sm_temp*noprime(Sm_RefBin*temp));
        std::complex<double> g2_cplx = g2_tensor.cplx();
        return g2_cplx;
    }
    
    std::complex<double> Mollow::calc_pop_bath(ITensor temp)
    {
        pop_bath_tensor = dag(noprime(Sm_RefBin*temp)) * noprime(Sm_RefBin*temp);            
        std::complex<double> pop_bath_cplx = pop_bath_tensor.cplx();
        return pop_bath_cplx;
    }

    double Mollow::norm(ITensor psi)
    {
        return (dag(psi)*psi).real();
    }

    void Mollow::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        std::string date = Mollow::get_date_as_string();
        
        std::ostringstream oss1, oss2, oss3;
        
        oss1 << "out/TLS_Mollow_timeevo" << date.c_str() << "_Om" << Omega << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_timeevo = oss1.str();
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "w");
        fprintf(datastream_timeevo,"#time \t population \t norm \n");
        fclose(datastream_timeevo);
        
        oss2 << "out/correlation_Mollow" << date.c_str() << "_Om" << Omega << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_correlation = oss2.str();
        datastream_correlation = fopen(file_identifier_correlation.c_str(), "w");
        fprintf(datastream_correlation, "#time \t Re[g1] \t Im[g1] \t Re[g2] \t Im[g2] \t pop_bath^2 \n");
        fclose(datastream_correlation);
        
        oss3 << "out/spectrum_Mollow" << date.c_str() << "_Om" << Omega << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbbins" << Nfb << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier_spectrum = oss3.str();
        datastream_spectrum = fopen(file_identifier_spectrum.c_str(), "w");
        fprintf(datastream_spectrum, "#omega \t Re[S] \t Im[S] \t bath_pop \n");
        fclose(datastream_spectrum);
        
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
    
    void Mollow::initialize_MPS()
    {
        // We need one tensor for each timestep and one for the system
        psi = MPS(ttotal+1);
        
        // system is at position A(Nfb+1) initially,timebins are initialized in vacuum
        ITensor MPSTensor = ITensor();

        for (int i=1;i<Nfb+1;i++)
        {
            MPSTensor = ITensor(bin[i]);
            MPSTensor.set(bin[i](1),1.0);
            psi.setA((i),MPSTensor);
        }
        
        MPSTensor = ITensor(bin[0]);
        MPSTensor.set(bin[0](1),init_11);
        MPSTensor.set(bin[0](2),init_22);
        printf("Initialize system with %.3f in ground and %.3f in upper level \n \n", init_11, init_22);
        psi.setA((Nfb+1),MPSTensor);
        
        for (int i=Nfb+2;i<ttotal+2;i++)
        {
            MPSTensor = ITensor(bin[i-1]);
            MPSTensor.set(bin[i-1](1),1.0);
            psi.setA((i),MPSTensor);
        }
        psi.orthogonalize();
    }
    
    void Mollow::set_up_MPO()
    {
        Hsys = ITensor(bin[0],prime(bin[0]));
        Hsys.set(bin[0](1), prime(bin[0](2)),-Cplx_i*Omega*dt);
        Hsys.set(bin[0](2), prime(bin[0](1)),-Cplx_i*Omega*dt);
        
        Hint = ITensor(bin[0],prime(bin[0]),bin[Nfb+1],prime(bin[Nfb+1]));
        Hfb = ITensor(bin[0],prime(bin[0]),bin[1],prime(bin[1]));

        for(int i=1;i<Nbin;i++)
        {            
            Hint.set(bin[0](2),prime(bin[0](1)),bin[Nfb+1](i),prime(bin[Nfb+1](i+1)), pow(dt,0.5)*Gamma*pow(i,0.5));
            Hint.set(bin[0](1),prime(bin[0](2)),bin[Nfb+1](i+1),prime(bin[Nfb+1](i)), (-1.)*pow(dt,0.5)*Gamma*pow(i,0.5));
            
            Hfb.set(bin[0](2),prime(bin[0](1)),bin[1](i),prime(bin[1](i+1)), (-1.0)*pow(dt,0.5)*Gamma*exp(-Cplx_i*phi)*pow(i,0.5));
            Hfb.set(bin[0](1),prime(bin[0](2)),bin[1](i+1),prime(bin[1](i)), pow(dt,0.5)*Gamma*exp(Cplx_i*phi)*pow(i,0.5));
        }
        
        ITensor identity_system = Mollow::identity(bin[0]);
        ITensor identity_int = Mollow::identity(bin[Nfb+1]);
        ITensor identity_fb = Mollow::identity(bin[1]);
        
        Hsys *= (identity_int*identity_fb);
        Hint *= identity_fb;
        Hfb *= identity_int;
        
        Hint += Hfb;
        
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
    
        Uevo = identity_system*identity_int*identity_fb + U_fb_1 + U_fb_2 + U_fb_3 + U_fb_4 + U_fb_5 + U_fb_6 + U_fb_7 + U_fb_8 + U_fb_9 + U_fb_10 ;
        
        // second order evolution
        // Hint += Hfb;
        // Uevo = identity_system*identity_int*identity_fb + Hsys + Hint + 0.5*mapprime(prime(Hint)*Hint,2,1);
    }
    
    void Mollow::timeevolve_system()
    {
        println("\n############################################################## \n \n");
        datastream_timeevo = fopen(file_identifier_timeevo.c_str(), "a");
        Mollow::setup_indexbin();
        Mollow::initialize_MPS();
        println("starting timeevolution of system... \n");
        printf("time %.3f of %.3f -- norm = %.10f -- occ = %.10f \n", 0.0*dt, ttotal, (1.0 - Mollow::norm(psi.A(1))), Mollow::occupation(psi.A(1), Site));
        Mollow::set_up_MPO();
        for(int j=Nfb; j<ttotal; j++)
        {
            Mollow::apply_Uevo(j);
            if(j<(ttotal-1))
            {
                // Uevo needs to act on future bin of MPS for next timestep:
                Uevo *= delta(bin[j+1],bin[j+2]);
                Uevo *= delta(prime(bin[j+1]),prime(bin[j+2]));
                Uevo *= delta(bin[j-Nfb+1],bin[j-Nfb+2]);
                Uevo *= delta(prime(bin[j-Nfb+1]),prime(bin[j-Nfb+2]));
            }
        }
        fclose(datastream_timeevo);
        println("\nfinished timeevolution of system.\n \n");
    }
    
    void Mollow::apply_Uevo(int timestep)
    {
        Mollow::swap_feedbackbin_to_system(timestep);
        
        ITensor U,S,V,temp;
        //shift ortho-center from feedbackbin to systembin to compute observable
        temp = psi.A(timestep)*psi.A(timestep+1);
        U = ITensor(bin[0],rightLinkInd(psi,timestep+1));
        svd(temp,U,S,V,args);
        psi.setA(timestep,V);
        psi.setA(timestep+1,U*S);
        if((timestep%save_every_nth_step)==0) 
        {
            fprintf(datastream_timeevo, "%f \t %e \t %e \t %e \t %e \n", (timestep-Nfb)*dt, Mollow::occupation(psi.A(timestep+1), Site), Mollow::ground_occupation(psi.A(timestep+1), Site), (Mollow::occupation(psi.A(timestep+1), Site)+Mollow::ground_occupation(psi.A(timestep+1), Site)), Mollow::norm(psi.A(timestep+1)));
        }
        
        if((timestep%print_every_nth_step)==0) 
        {
            printf("time %.3f of %.3f -- norm = %e -- occ = %e \n",(timestep-Nfb)*dt, (ttotal-Nfb)*dt, (1.0 - Mollow::norm(psi.A(timestep+1))), Mollow::occupation(psi.A(timestep+1), Site));
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
        
        Mollow::swap_feedbackbin_to_origin(timestep);
    }
    
    void Mollow::swap_feedbackbin_to_system(int timestep)
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
    
    void Mollow::swap_feedbackbin_to_origin(int timestep)
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
    
    void Mollow::calc_correlation()
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
        Sm_RefBin = Mollow::annihilator(bin[referencetimebin]);
        
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
            g1_bath[k] = Mollow::calc_g1_bath(s, temp);
            g2_bath[k] = Mollow::calc_g2_bath(s, temp);
            pop_bath[k] = Mollow::calc_pop_bath(temp);
            
            if((i%save_every_nth_step)==0) 
            {
                fprintf(datastream_correlation,"%e \t %e \t %e \t %e \t %e \t %e \t %e \n",dt*k,g1_bath[k].real(),g1_bath[k].imag(), g2_bath[k].real(), g2_bath[k].imag(), pop_bath[k].real(), pop_bath[k].real()*pop_bath[k].real());
            }
            
            if((i%print_every_nth_step)==0) 
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
    
    void Mollow::calc_spectrum()
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


        // timeloop endet bei i=ttotal-1. Dann swape ich zurück und der Timebin ist bei 
        // psi.A(ttotal-1-Nfb+1) = psi.A(ttotal-Nfb).
        // Aber dort ist das orthocenter nicht: das ist dann bei psi.A(ttotal-Nfb+1).
        // dh. ich muss erstmal das orthocenter auf psi.A(ttotal-Nfb) schieben.

