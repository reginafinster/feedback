#include "feedback.h"


namespace itensor{
    
    Feedback::Feedback(std::string infile)
    {
        InputGroup input = InputGroup(infile,"input");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        dt = input.getReal("time_step");
        t_end = input.getReal("time_end");
        Omega = input.getReal("Omega_TLS");
        Gamma = input.getReal("Gamma");
        phi = input.getReal("feedbackphase_in_units_of_phi",0);
        Nbin = input.getInt("Dim_of_bath");
        Nfb = input.getInt("feedback_time");
        init_cond = input.getString("initial_condition");
        ttotal = t_end;
        ttotal += Nfb;
        args = Args("Cutoff=",cutoff,"Maxm=",maxm);
        phi *= 3.141592654;
        Feedback::setup_dir_and_filename();
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

    double Feedback::norm(ITensor psi, IndexType type)
    {
        return (dag(psi)*psi).real();
    }

    void Feedback::setup_dir_and_filename()
    {
        mkdir("./out", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        
        std::string date = Feedback::get_date_as_string();
        std::ostringstream oss;
        oss << "out/TLS_Feedback_" << date.c_str() << "_Om" << Omega << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_fbtime" << (Nfb*dt) << "_phase" << (phi/3.141592654) << ".dat";
        file_identifier = oss.str();
        datastream = fopen(file_identifier.c_str(), "w");
        fprintf(datastream,"#time \tS+S- \t norm \n");
        fclose(datastream);
    }
    
    void Feedback::setup_indexbin()
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
    
    void Feedback::initialize_MPS()
    {
        // We need one tensor for each timestep and one for the system
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
            MPSTensor.set(bin[0](2),1.0);
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
    
    void Feedback::set_up_MPO_with_feedback()
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
        
        ITensor identity_system = Feedback::identity(bin[0]);
        ITensor identity_int = Feedback::identity(bin[Nfb+1]);
        ITensor identity_fb = Feedback::identity(bin[1]);
        
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
    
    void Feedback::timeevolve_system_with_feedback()
    {
        datastream = fopen(file_identifier.c_str(), "a");
        
        Feedback::setup_indexbin();
        Feedback::initialize_MPS();
        printf("time %.3f of %.3f -- norm=%.10f -- occ=%.10f \n", 0.0*dt, t_end, (1.0 - Feedback::norm(psi.A(1), Site)), Feedback::TLS_occupation(psi.A(1), Site));
        Feedback::set_up_MPO_with_feedback();
        for(int j=Nfb; j<ttotal; j++)
        {
            Feedback::apply_Uevo_with_feedback(j);
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
    
    void Feedback::apply_Uevo_with_feedback(int timestep)
    {
        Feedback::swap_feedbackbin_to_system(timestep);
        
        ITensor U,S,V,temp;
        //shift ortho-center from feedbackbin to systembin to compute observable
        temp = psi.A(timestep)*psi.A(timestep+1);
        U = ITensor(bin[0],rightLinkInd(psi,timestep+1));
        svd(temp,U,S,V,args);
        psi.setA(timestep,V);
        psi.setA(timestep+1,U*S);
        fprintf(datastream, "%f \t %e \t %e \n", (timestep-Nfb)*dt, Feedback::TLS_occupation(psi.A(timestep+1), Site), (1.0 - Feedback::norm(psi.A(timestep+1), Site)));
        if ((timestep%10)==0) 
        {
            printf("time %.3f of %.3f -- norm=%.10f \n",(timestep-Nfb)*dt,t_end*1., (1.0 - Feedback::norm(psi.A(timestep+1), Site)));
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
        
        Feedback::swap_feedbackbin_to_origin(timestep);
    }
    
    void Feedback::swap_feedbackbin_to_system(int timestep)
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
    
    void Feedback::swap_feedbackbin_to_origin(int timestep)
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

}//eof namespace




