#include "spinchain_feedback.h"


namespace itensor{
    
    Feedback::Feedback(std::string infile)
    {
        InputGroup input = InputGroup(infile,"input");
        cutoff = input.getReal("svdcutoff");
        maxm = input.getInt("maxnumberofSV");
        steadybins = input.getInt("steadybins",10);
        dt = input.getReal("time_step");
        t_end = input.getReal("time_end");
        N = input.getInt("N",4);
        J = input.getReal("J");
        Jz = input.getReal("Jz");
        Omega = input.getReal("Omega_TLS");
        Gamma = input.getReal("Gamma");
        phi_start = input.getReal("phi_start",-0.1);
        phi_stepsize = input.getReal("phi_stepsize",0.1);
		phi_steps = input.getInt("phi_steps",13);
        Nbin = input.getInt("Dim_of_bath");
        Nfb = input.getInt("feedback_time");
        enable_feedback = input.getInt("enable_feedback",0);
        read_in_shifts = input.getInt("read_in_shifts",0);
        initial = input.getInt("initial",0);
        init_pos = input.getInt("init_pos",1);
        init_pos1 = input.getInt("init_pos1",1);
        init_pos2 = input.getInt("init_pos2",1);
        init_pos3 = input.getInt("init_pos3",1);
        init_pos4 = input.getInt("init_pos4",1);
        init_pos5 = input.getInt("init_pos5",1);
        neel_state_first_up = 0;
        neel_state_first_down = 0;
        allup_state = 0;
        one_excitation = 0;
        various_excitations = 0;
        if(initial == 1) allup_state = 1;
		if(initial == 2) neel_state_first_up = 1;
        if(initial == 3) neel_state_first_down = 1;
        if(initial == 4) one_excitation = 1;
        if(initial == 5) various_excitations = 1;
        ttotal = t_end/dt;
        ttotal += Nfb;
        args = Args("Cutoff=",cutoff,"Maxm=",maxm);
    }
    
    std::string Feedback::get_date_as_string()
    {
        time_t t = time(0); 
        struct tm *now = localtime(&t);
        std::stringstream datestring;
        datestring << now->tm_mday << '-' << (now->tm_mon + 1) << '-' << (now->tm_year + 1900);
        return datestring.str();
    }
    
    ITensor Feedback::identity(Index& s)
    {
        ITensor my_identity = ITensor(s,prime(s));
        for(int i=1; i<=s.m(); i++)
        {
            my_identity.set(s(i),prime(s)(i),1.0);
        }
        return my_identity;
    }
    
    ITensor Feedback::MySz(Index& s)
	//this function creates the operator on the site i of the given Index
	//input: SiteSet and number i (type int) of site
	//return value: Operator Sz on site i (type ITensor)
	{
		ITensor MySz = ITensor(s,prime(s));
		MySz.set(s(1),prime(s(1)),0.5);
		MySz.set(s(2),prime(s(2)),-0.5);
		return MySz;
	}

	ITensor Feedback::MySp(Index& s)
	{
		ITensor MySp = ITensor(s,prime(s));
		MySp.set(s(2),prime(s(1)),1.0);
		return MySp;
	}

	ITensor Feedback::MySm(Index& s)
	{
		ITensor MySm = ITensor(s,prime(s));
		MySm.set(s(1),prime(s(2)),1.0);
		return MySm;
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
    
    double Feedback::occupation(MPS& MPS_chain, int i)
	//this function calculates the magnetization of an antiferromagnetic spin chain.
	//z-value of spin is counted positive for spin up and negative for spin down
	{
		MPS_chain.position(i);
		ITensor SpSm = ITensor(s[i], prime(s[i]));
        SpSm.set(s[i](2), prime(s[i](2)), 1.0);
		double occ = real((dag(prime(MPS_chain.A(i),s[i],1)) * SpSm * MPS_chain.A(i)).cplx());
		return occ;
	}
    
    double Feedback::stag_magn(MPS& MPS_chain)
	//this function calculates the magnetization of an antiferromagnetic spin chain.
	//z-value of spin is counted positive for spin up and negative for spin down
	{
		int N = MPS_chain.N();
		double stagMagn = 0.0;     
		
		//loop through chain and add single magnetization of each site
		for (int i=1;i<=N;i++)
		{
			MPS_chain.position(i);
			Index temp = s[i]; //TODO evtl ausserhalb des loops deklarieren? auch mag_on...??
			Complex mag_one_site = (dag(prime(MPS_chain.A(i),s[i],1)) * MySz(temp) * MPS_chain.A(i)).cplx();
			if (i%2 == 0) {stagMagn += real(mag_one_site)/N;}
			else if (i%2 == 1) {stagMagn -= real(mag_one_site)/N;}
		}    
		return stagMagn;
	}
	
	double Feedback::integrate_detector(MPS& MPS_bins,int startbin,int lastbin)
	{
		double erg = 0.0;
		for(int i=startbin;i<lastbin+1;i++)
		{
			MPS_bins.position(i);
			ITensor SpSm = ITensor(bin[i], prime(bin[i]));
			SpSm.set(bin[i](2), prime(bin[i](2)), 1.0);
			erg += real((dag(prime(MPS_bins.A(i),bin[i],1)) * SpSm * MPS_bins.A(i)).cplx());
		}
		return erg;
	}	
	
	double Feedback::transient_detector(ITensor MPSTensor)
	{
		Index binindex = findtype(MPSTensor, Site);
		ITensor SpSm = ITensor(binindex, prime(binindex));
		SpSm.set(binindex(2), prime(binindex(2)), 1.0);
		double erg = real((dag(prime(MPSTensor,binindex,1)) * SpSm * MPSTensor).cplx());
		return erg;
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
        std::ostringstream oss1;
        std::ostringstream oss2;
        std::ostringstream oss3;
        oss << "out/Spinchain_Feedback_" << date.c_str() << "_N" << N << "_J" << J << "_Gam"<< Gamma << "_dt" << dt << "_cutoff" << cutoff << "_maxm" << maxm << "_init" << initial << "_pos" << init_pos;
		if(read_in_shifts) oss << "_shifted";
        if(enable_feedback) oss << "_Nfb" << Nfb << "_fbp" << inphi;
        oss1 << oss.str();
        oss2 << oss.str();
        oss3 << oss.str();
        //std::cout << oss1.str(); 
        
		oss << "_occupation.dat";
        file_identifier = oss.str();
        datastream = fopen(file_identifier.c_str(), "w");
        fprintf(datastream,"#time");
        for (int i=N;i>0;i--) fprintf(datastream,"\t occupation %i",i);//standard xmgrace plot should give last site (e.g. for TLS BM etc)
        fprintf(datastream," \t norm \n");
        fclose(datastream);
        
        oss1 << "_stagmag.dat";
        file_identifier1 = oss1.str();
        stagmagstream = fopen(file_identifier1.c_str(), "w");
        fprintf(stagmagstream,"#time \t stagmag \n");
        fclose(stagmagstream);
        
        oss2 << "_steady.dat";
        file_identifier2 = oss2.str();
        steadystream = fopen(file_identifier2.c_str(), "w");
        fprintf(steadystream,"# J \t\t Gamma \t\t tau \t\t phi \t\t initial \t\t init_pos \t\t shifted \t\t Occupation %i", N);
        for(int i = N-1;i>0;i--) fprintf(steadystream,"\t\t Occupation %i", i);
        fprintf(steadystream,"\t\t detector \t\t trapped_bins \t\t detector+all_occ \t\t steadycheck \t\t norm \n");
        fclose(steadystream);
                
        oss3 << "_detector.dat";
        file_identifier3 = oss3.str();
        detectorstream = fopen(file_identifier3.c_str(), "w");
        fprintf(detectorstream,"#time \t detector \n");
        fclose(detectorstream);
    }
    
    void Feedback::setup_indices()
    {
        // TODO first entry of vector is the system - not anymore! wie muss der initialisiert werden? trotzdem so oder?
        
        bin = std::vector<Index>(ttotal+2);
        for (int i=1;i<(ttotal+2);i++)
        {
            bin.at(i) = Index(nameint("k",i),Nbin,Site);
        }
        
        s = std::vector<Index>(N+1); //Index for Sites (Spins)
		for (int i=1; i<N+1; i++)
		{
			s.at(i) = Index(nameint("s",i), 2, Site);
		}        
    }
    
    void Feedback::setup_shifts(double* h)
    {
		std::string line;
		std::ifstream shiftsfile ("shifts.dat");
		if (shiftsfile.is_open())
		{
			int m=0;
			while ( getline (shiftsfile,line) )
			{
				if((m>0) && (m<N+1)) //skip first line to match RK Benchmark
				{
					h[m]=stod(line);
					printf("h%i = %f\n",m,h[m]);
				}
				m++;
			}
		}else 
		{
			printf("Unable to open shiftsfile\n");
		}
		shiftsfile.close();
    }
    
    void Feedback::initialize_MPS()
    {	
		println("\n set up two MPS\n");
		ITensor MPSTensor = ITensor();
        // system is initialized
		MPS_chain = MPS(N);
        if(one_excitation){
			for (int i=1;i<N+1;i++)
			{
				MPSTensor = ITensor(s[i]);
				if(i==init_pos) {MPSTensor.set(s[i](2),1.0);}
				else {MPSTensor.set(s[i](1),1.0);}
				MPS_chain.setA((i),MPSTensor);
			}
		}else if(various_excitations)
        {
            for (int i=1;i<N+1;i++)
			{
				MPSTensor = ITensor(s[i]);
				if(i==init_pos)  {MPSTensor.set(s[i](2),1.0);}
				else if(i==init_pos1) {MPSTensor.set(s[i](2),1.0);}
				else if(i==init_pos2) {MPSTensor.set(s[i](2),1.0);}
				else if(i==init_pos3) {MPSTensor.set(s[i](2),1.0);}
				else if(i==init_pos4) {MPSTensor.set(s[i](2),1.0);}
				else if(i==init_pos5) {MPSTensor.set(s[i](2),1.0);}
				else {MPSTensor.set(s[i](1),1.0);}
				MPS_chain.setA((i),MPSTensor);
			}
        }
        else
		{
			int even = 1;
			int odd = 1;
			if(neel_state_first_up) odd = 2;
			if(neel_state_first_down) even = 2;
			if(allup_state) {even = 2; odd = 2;}
			
			for (int i=1;i<N+1;i++)
			{
				MPSTensor = ITensor(s[i]);
				if(i%2==0) {MPSTensor.set(s[i](even),1.0);}
				else {MPSTensor.set(s[i](odd),1.0);}
				MPS_chain.setA((i),MPSTensor);
			}
		}
        PrintData(MPS_chain);
        MPS_chain.orthogonalize();

        // We need one tensor for each timestep and one for the system
        printf("ttotal = %i\n",ttotal);
        MPS_bins = MPS(ttotal+2);
        // timebins are initialized in vacuum
        // last site of system is at position A(Nfb+1) initially
        MPS_bins.setA(1,MPS_chain.A(N));
		for (int i=1;i<ttotal+2;i++)
		{
			MPSTensor = ITensor(bin[i]);
			MPSTensor.set(bin[i](1),1.0);
			MPS_bins.setA((i+1),MPSTensor);
		}
		MPS_bins.orthogonalize();
		//Print(MPS_bins);
        if(enable_feedback)
        {      
			ITensor U,S,V,temp;
			Index systemtimebin = findtype(MPS_bins.A(1),Site);
			Index systemlink = commonIndex(MPS_chain.A(N-1),MPS_chain.A(N),Link);
			for(int i=1; i<ttotal+1; i++)
			{
				temp = MPS_bins.A(i)*MPS_bins.A(i+1);
				if(i<Nfb+1)
				{
					U = ITensor(systemtimebin, systemlink, commonIndex(MPS_bins.A(i+1),MPS_bins.A(i+2),Link));
				}else
				{
					U = ITensor(bin[i],commonIndex(MPS_bins.A(i+1),MPS_bins.A(i+2),Link));
				}
				svd(temp,U,S,V,args);
				MPS_bins.setA(i,V);
				MPS_bins.setA(i+1,U*S);
				//Print(MPS_bins);
			}
		}
        //Print(MPS_bins);
	//in any case the orthocenter should now be on the bin containing s_N.
    }
    
    void Feedback::set_up_MPO()
    {
		ITensor identity_sys_last = Feedback::identity(s[N]);
		ITensor identity_sys_before_last = Feedback::identity(s[N-1]);
		ITensor identity_bin = Feedback::identity(bin[1]);
		ITensor identity_Nfb = Feedback::identity(bin[Nfb+1]);
		ITensor L,R,identityN;

		if(enable_feedback)
		{
			println("\n set up MPO with feedback\n");

			//Hsys = ITensor(bin[0],prime(bin[0])); //treibung erstmal nicht vorgesehen. später vielleicht?
			//Hsys.set(bin[0](1), prime(bin[0](2)),-Cplx_i*Omega*dt);
			//Hsys.set(bin[0](2), prime(bin[0](1)),-Cplx_i*Omega*dt);
			
			//Feedback-Part
			Hint = ITensor(s[N],prime(s[N]),bin[Nfb+1],prime(bin[Nfb+1]));//bin[0] is now s[N] of other MPS
			Hfb = ITensor(s[N],prime(s[N]),bin[1],prime(bin[1]));

			for(int i=1;i<Nbin;i++)
			{            
				Hint.set(s[N](2),prime(s[N](1)),bin[Nfb+1](i),prime(bin[Nfb+1](i+1)), pow(dt,0.5)*Gamma*pow(i,0.5));//Kommunikation mit "zukunft"
				Hint.set(s[N](1),prime(s[N](2)),bin[Nfb+1](i+1),prime(bin[Nfb+1](i)), (-1.)*pow(dt,0.5)*Gamma*pow(i,0.5));        
				Hfb.set(s[N](2),prime(s[N](1)),bin[1](i),prime(bin[1](i+1)), (-1.0)*pow(dt,0.5)*Gamma*exp(-Cplx_i*phi)*pow(i,0.5));//Kommunikation mit fb bin
				Hfb.set(s[N](1),prime(s[N](2)),bin[1](i+1),prime(bin[1](i)), pow(dt,0.5)*Gamma*exp(Cplx_i*phi)*pow(i,0.5));
			}
			
			//Hsys *= (identity_bin*identity_Nfb);
			Hint *= identity_bin;
			Hfb *= identity_Nfb;
			
			Hint += Hfb;

		}else
		{
			println("\n set up MPO without feedback\n");
			Hint = ITensor(s[N],prime(s[N]),bin[1],prime(bin[1]));//UNKLAR muss das angepasst werden? Nein, weil es hier keinen FB-Teil der Timebin-chain gibt.
			
			for (int i=1; i<Nbin; i++)
			{
				Hint.set(s[N](2),prime(s[N](1)),bin[1](i),prime(bin[1](i+1)), Gamma*pow(2.*dt,0.5)*pow(i,0.5));//factor sqrt(2), since only one output, not 2 like in the fb case.
				Hint.set(s[N](1),prime(s[N](2)),bin[1](i+1),prime(bin[1](i)), (-1.0)*Gamma*pow(2.*dt,0.5)*pow(i,0.5));
			}
		}

		Uevo = Hint*identity_sys_before_last; //expand to full indexspace for summation later

		//Chain-MPO
        MPO H_E = MPO(N);
		MPO H_O = MPO(N);
		double h[N];
		if(read_in_shifts)
		{
			Feedback::setup_shifts(h);
		}
		for(int i=1;i<N;i++)
		{
			//define Heisenberg-Hamiltonian
			ITensor H_temp  = 4.*Jz*Feedback::MySz(s[i])*Feedback::MySz(s[i+1]);
					H_temp += 4.*0.5*J*Feedback::MySp(s[i])*Feedback::MySm(s[i+1]);
					H_temp += 4.*0.5*J*Feedback::MySm(s[i])*Feedback::MySp(s[i+1]);
					if(read_in_shifts)
					{
						H_temp += h[i]*Feedback::MySz(s[i])*Feedback::identity(s[i+1]);
						H_temp += h[i+1]*Feedback::MySz(s[i+1])*Feedback::identity(s[i]); 
						if(i==1) H_temp += h[1]*Feedback::MySz(s[1])*Feedback::identity(s[2]); //otherwise, rim-shifts have only half the effect of middle shifts
						if(i==N-1) H_temp += h[N]*Feedback::MySz(s[N])*Feedback::identity(s[N-1]);
					}
			//H_E and H_0 each get the full set of indices s, but now they are factorized into MPOs
			if (i == N-1) //last site is special case
			{	
				H_temp *= dt;
				//PrintData(H_temp);
				if(enable_feedback)
				{
					H_temp *= identity_bin*identity_Nfb;
				}else
				{
					H_temp *= identity_bin;
				}
				H_temp += Uevo;
				//PrintData(H_temp);
				Cplx cplx_dt = (-1.0)*Cplx_i;
				H_temp = expHermitian(H_temp,cplx_dt);//UNKLAR mögl Fehlerquelle: da H hermitesch ist kann das bleiben oder? oder vllt besser auch mit Taylor? 
				
				if(enable_feedback){
					L = ITensor(s[i],prime(s[i]));
					R = ITensor(s[i+1],prime(s[i+1]),bin[1],prime(bin[1]),bin[Nfb+1],prime(bin[Nfb+1])); //all timebin-indices on last site
				}else
				{
					L = ITensor(s[i],prime(s[i]));
					R = ITensor(s[i+1],prime(s[i+1]),bin[1],prime(bin[1])); 
				}
				factor(H_temp,L,R);
				
				//PrintData(L);				
				//PrintData(R);				
				
				if (i%2 == 0)
				{
					H_E.setA(i,L);
					H_E.setA(i+1,R);
				}
				else if (i%2 ==1)
				{
					H_O.setA(i,L);
					H_O.setA(i+1,R);
				}
			}else
			{
				H_temp *= dt;
				Cplx cplx_dt = (-1.0)*Cplx_i;
				H_temp = expHermitian(H_temp,cplx_dt); //UNKLAR mögl Fehlerquelle: da H hermitesch ist kann das bleiben oder? oder vllt besser auch mit Taylor? 
			
				ITensor L(s[i],prime(s[i])), R(s[i+1],prime(s[i+1]));
				factor(H_temp,L,R); //UNKLAR warum hier nicht svd??
				
				//PrintData(L);				
				//PrintData(R);				
				
				if (i%2 == 0)
				{
					H_E.setA(i,L);
					H_E.setA(i+1,R);
				}
				else if (i%2 == 1)
				{
					 H_O.setA(i,L);
					 H_O.setA(i+1,R);
				}
			}
		printf("%i\n",i);
		//PrintData(H_O);
		//PrintData(H_E);
		}
		
		//the first Tensor of H_E is always empty and last of H_E (N is even) or H_O (N is odd)
		//is still empty and needs to be filled with Identity
		ITensor identity1 = Feedback::identity(s[1]);//TODO anpassen
		H_E.setA(1,identity1);
		
		if(enable_feedback)
		{
			identityN = Feedback::identity(s[N])*Feedback::identity(bin[1])*Feedback::identity(bin[Nfb+1]);//last site needs all indices of the respective other, so system and timebin-indices
		}else
		{
			identityN = Feedback::identity(s[N])*Feedback::identity(bin[1]);//last site needs all indices of the respective other, so system and timebin-indices
		}
		//PrintData(identityN);
		if (N%2==0)
		{
			H_E.setA(N,identityN);
		}
		else if (N%2==1)
		{
			H_O.setA(N,identityN);
		}
		
		//PrintData(H_O);
		//PrintData(H_E);
				
		//multiply H_E and H_O, so that tensors on the same index are contracted
		//U is now the operator for one Trotter-step
		MPO_chain = Feedback::multiply_two_nonherm_MPO(H_E,H_O); //TODO Funktion einfügen
		Print(MPO_chain);
    }
    
    MPS Feedback::apply_nonherm_MPO(MPO& Uevo, MPS& psi, MPS& res, Args args)
	// Beachte bei der Eingabe: MPS und MPO müssen mit den selben Site-Indizes versehen sein. Bei beiden muss der 
	// Zeilenindex ungeprimed und der Spaltenindex geprimed sein. 
	// Ist dem nicht so, wird faelschlicherweise mit dem transponierten MPO multipliziert!
	{
		// set up empty MPS for result
		int N = psi.N();
		Uevo.position(1);
		psi.position(1);

		// prime site indices of MPS in order to achieve correct multiplication with non-hermitian MPO
		psi.mapprime(0,1,Site);
		
		// declare indices and tensors
		Index si_left;
		Index si_right;
		Index li_psi;
		Index li_Uevo;
		Index li_res;
		ITensor *temp = new ITensor;
		ITensor *U = new ITensor;
		ITensor *S = new ITensor;
		ITensor *V = new ITensor;
		
		 // loop through sites. First and last site are special cases
		for (int i = 1; i<N; i++)
		{
			Index si_left = mapprime(commonIndex(psi.A(i), Uevo.A(i), Site),1,0);
			Index si_right = mapprime(commonIndex(psi.A(i+1), Uevo.A(i+1), Site),1,0);

			// contract tensors depending on their position in the chain
			if (i == 1)
			{
				// define additional link indices
				li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
				li_Uevo = commonIndex(Uevo.A(i+1), Uevo.A(i+2), Link);
				// contract Tensors and decompose them again
				*temp = psi.A(i) * Uevo.A(i) * psi.A(i+1) * Uevo.A(i+1);
				*U = ITensor(si_left);
				*V = ITensor(si_right, li_psi, li_Uevo);
			}
			else if (i == (N-1))
			{
				// define additional link indices
				li_res = commonIndex(res.A(i-1), res.A(i), Link);
				// contract Tensors and decompose them again
				*temp = res.A(i) * psi.A(i+1) * Uevo.A(i+1);
				*U = ITensor(si_left, li_res);
				*V = ITensor(si_right);
			}
			else
			{
				// define additional link indices
				li_res = commonIndex(res.A(i-1), res.A(i), Link);
				li_psi = commonIndex(psi.A(i+1), psi.A(i+2), Link);
				li_Uevo = commonIndex(Uevo.A(i+1), Uevo.A(i+2), Link);
				// contract Tensors and decompose them again
				*temp = res.A(i) * psi.A(i+1) * Uevo.A(i+1);
				*U = ITensor(si_left, li_res);
				*V = ITensor(si_right, li_psi, li_Uevo);
			}
			svd(*temp, *U, *S, *V, args);
			*V = (*V) * (*S);
			//factor(temp, U, V, args);
			// store Tensors into resulting MPS:
			res.setA(i,*U);
			res.setA(i+1,*V);
		}//EOF site loop
		
		// restore original state of psi
		psi.mapprime(1,0,Site); //UNKLAR: psi wird doch sowieso von res überschrieben in main oder??!?
		delete U;
		delete V;
		delete S;
		delete temp;
		
		return res;
	}
	
	MPO Feedback::multiply_two_nonherm_MPO(const MPO& MPO_A, const MPO& MPO_B)
	// Beachte bei der Eingabe: Die beiden MPOs müssen mit den selben Site-Indizes versehen sein. Bei beiden muss der 
	// Zeilenindex ungeprimed und der Spaltenindex geprimed sein. 
	{
		// set up empty MPS for result
		int N = MPO_A.N();
		int K = MPO_B.N(); 
		if (N!=K) 
		{ 
			std::cout << "In function multiply_two_nonherm_MPO: size of MPOs does not fit. Abort" << std::endl; 
		}
		
		MPO A = MPO_A;
		MPO B = MPO_B;
		MPO res = MPO(N);
		//PrintData(A);
		//PrintData(B);
		
		A.position(1);
		B.position(1);
		
		B.primeall();
		
		Index si_A;
		Index li_res;
		ITensor temp, fork, U, S, V;
		auto IndexTypeSiteNoPrime = [](Index const& i) { return (i.type() == Site && i.primeLevel() == 0); };
		
		// loop through sites. First and last site are special cases
		for (int i = 1; i<=N; i++)
		{
			if (i == 1)
			{ 
				si_A = findindex(A.A(i), IndexTypeSiteNoPrime);
				temp = A.A(i) * B.A(i);
				U = ITensor(si_A, prime(si_A,2));
				svd(temp, U, S, V);
				U = U * S;
				res.setA(i,U);
			}
			else if (i == N)
			{
				si_A = findindex(A.A(i), IndexTypeSiteNoPrime);
				li_res = commonIndex(U, V, Link);
				temp = V * A.A(i) * B.A(i);
				res.setA(i,temp);
			}
			else
			{
				si_A = findindex(A.A(i), IndexTypeSiteNoPrime);
				li_res = commonIndex(U, V, Link);
				temp = V * A.A(i) * B.A(i);
				U = ITensor(si_A, prime(si_A,2), li_res);
				svd(temp, U, S, V);
				U = U * S;
				res.setA(i,U);
			}
			
		} // EOF site loop
		
		res.mapprime(2,1,Site);
		res.orthogonalize();
		return res;  
	}
    
    void Feedback::timeevolve_system()
    // past=j as MPS starts at 1 and bin[] at 0
    // system is always at psi(j+1), future at psi(j+2)
    // Uevo acts on system and future timebin
    // thus we contract these two tensors, apply Uevo and decompose
    // with U(futurebin, link to past), V(systembin, link to future)
    // and set orthocenter to V 
    // U becomes new past bin in MPS
    // V becomes new systembin at next timestep in MPS
    {
		for(int k = 0; k<phi_steps;k++)
		{
			phi = phi_start + k*phi_stepsize;
			inphi = phi;
			printf("start of interation for phi = %f pi\n\n",inphi);
			phi *= 3.141592654;
			Feedback::setup_dir_and_filename();
			
			datastream = fopen(file_identifier.c_str(), "a");
			stagmagstream = fopen(file_identifier1.c_str(), "a");
			detectorstream = fopen(file_identifier3.c_str(), "a");
			ITensor temp;
			trans_det = 0.0;
            dissipated_signal = 0.0;
			Feedback::setup_indices();
			Feedback::initialize_MPS();
			printf("MPS is initialized.\n");
			printf("time %.3f of %.3f -- norm_last_site=%.10f -- occ_last_site=%.10f \n", 0.0*dt, t_end, (1.0 - Feedback::norm(MPS_chain.A(N), Site)), Feedback::occupation(MPS_chain, N));//TODO was wollen wir plotten? stagmag?
			
			Feedback::set_up_MPO();
			int lastbin = 0;
			if(enable_feedback)
			{
				printf("start timeevolution of system with feedback\n");
				for(int j=Nfb; j<ttotal; j++)
				//for(int j=Nfb; j<Nfb+3; j++)
				{
					Feedback::apply_Uevo_with_feedback(j);
					//print("applied.\n");
					if(j<(ttotal-1))
					{
						// Uevo needs to act on future bin of MPS for next timestep:
						temp = MPO_chain.A(N);
						//Print(temp);
						temp = temp*delta(bin[j+1],bin[j+2])*delta(prime(bin[j+1]),prime(bin[j+2]))*delta(bin[j-Nfb+1],bin[j-Nfb+2])*delta(prime(bin[j-Nfb+1]),prime(bin[j-Nfb+2]));
						//Print(temp);
						MPO_chain.setA(N,temp);
					}
				lastbin = j;
				}
			}else
			{
				printf("start timeevolution of system without feedback\n");
				for(int j=1; j<ttotal; j++)
				//for(int j = 1; j<5; j++)//iteration starts with 2, since 1 and two are special cases, done in one step.
				{
					Feedback::apply_Uevo_without_feedback(j);
					//if(j < 3 ) PrintData(MPS_chain);
					if(j<(ttotal-1))
					{
						temp = MPO_chain.A(N);
						//Print(temp);
						temp = temp*delta(bin[j],bin[j+1])*delta(prime(bin[j]),prime(bin[j+1]));
						MPO_chain.setA(N,temp);
						//Print(MPO_chain);
					}
				lastbin = j;
				}
			}
			printf("finished timevolution\n");
			steadystream = fopen(file_identifier2.c_str(), "a");
			fprintf(steadystream, "%f \t\t %f \t\t %f \t\t %f \t\t %i \t\t %i \t\t %i",J,Gamma,dt*Nfb,inphi,initial,init_pos,read_in_shifts);
			for (int i = N; i>0; i--)
			{
			MPS_chain.position(i);
			fprintf(steadystream, "\t\t %.10e ",Feedback::occupation(MPS_chain, i));
			}
			double steadycheck = Feedback::integrate_detector(MPS_bins,lastbin-Nfb-steadybins,lastbin-Nfb);
			double detector_int = Feedback::integrate_detector(MPS_bins,1,lastbin-Nfb);
			double trapped_bins = Feedback::integrate_detector(MPS_bins,1,lastbin)-detector_int;
			double controlsum = detector_int+trapped_bins;
			for(int i=1;i<N+1;i++) controlsum += Feedback::occupation(MPS_chain, i);
			fprintf(steadystream,"\t\t %f \t\t %f \t\t %f \t\t %f \t\t %f \n",detector_int,trapped_bins,controlsum,steadycheck,Feedback::norm(MPS_chain.A(1), Site));
			fclose(steadystream);
			fclose(datastream);
			fclose(stagmagstream);
			fclose(detectorstream);
			printf("end of iteration for phi = %f pi\n\n",inphi);
		}
    }
    
    void Feedback::apply_Uevo_with_feedback(int timestep)
    {
		//int i = 1;
		//printf("hier%i\n",i);i++;
        Feedback::swap_feedbackbin_to_system(timestep);
        
        ITensor U,S,V,temp;
        MPS_bins.position(timestep);
        Index currenttimebin = findtype(MPS_bins.A(timestep+2),Site);
        //printf("hier%i\n",i);i++;
        
        //contract system-bin with interaction and feedback bin
        temp = noprime(MPS_bins.A(timestep)*MPS_bins.A(timestep+1)*MPS_bins.A(timestep+2));//TODO hier muss dieses triplet in die MPS_chain geschrieben werden und dann der MPO angewendet werden.
        //write this into the spin chain
        //Print(MPS_chain);
        MPS_chain.setA(N,temp);
		//Print(MPS_chain);
		//printf("hier%i\n",i);i++;
        
        //apply the MPO
        MPS *res = new MPS;
		*res = MPS(N);
        MPS_chain = Feedback::apply_nonherm_MPO(MPO_chain, MPS_chain, *res, args);
        delete res;
        if(timestep % 50 == 0) printf("time %.3f of %.3f -- norm_first_site=%.10f -- stagmag=%.10f -- occ_last_site=%.10f -- detector=%.10f\n",
										(timestep-Nfb)*dt, t_end, Feedback::norm(MPS_chain.A(1), Site), Feedback::stag_magn(MPS_chain), Feedback::occupation(MPS_chain, N),trans_det);
        if(timestep % 1 == 0)
        {
            fprintf(datastream, "%.3f ", (timestep-Nfb+1)*dt);
            for (int i=N;i>0;i--) fprintf(datastream, "\t %.5e ", Feedback::occupation(MPS_chain, i));
            fprintf(datastream, "\t %.10f \n", Feedback::norm(MPS_chain.A(1), Site));
            fprintf(stagmagstream, "%.3f \t %.5e \n", (timestep-Nfb+1)*dt, Feedback::stag_magn(MPS_chain));
            fprintf(detectorstream, "%.3f \t %.5e \t %.5e \n", (timestep-Nfb+1)*dt, trans_det, dissipated_signal);
        }
        //printf("hier%i\n",i);i++;
        
        //decompose and write the triplet back into the bin-chain (different positions)
        temp = MPS_chain.A(N);
        // bring system to A(timestep+2), thus to the left
        U = ITensor(s[N],leftLinkInd(MPS_bins,timestep+3),leftLinkInd(MPS_chain,N));
        svd(temp,U,S,V,args);
        V *= S;
        //Print(V);
        //printf("hier%i\n",i);i++;
        MPS_bins.setA(timestep+2,U);
            
        // bring current timebin one position to the left (which is middle position)
        // and feedbacktimebin back to its position at A(timestep) 
        U = ITensor(currenttimebin,commonIndex(U,V,Link));
        temp=V;
        svd(temp,U,S,V,args);
        MPS_bins.setA(timestep+1,U);
        MPS_bins.setA(timestep,V*S);
        dissipated_signal = Feedback::transient_detector(V*S);
        trans_det += Feedback::transient_detector(V*S);
        //printf("hier%i\n",i);i++;
        //Print(MPS_bins);
        Feedback::swap_feedbackbin_to_origin(timestep);
    }
    
    void Feedback::apply_Uevo_without_feedback(int timestep)
    {	
		int i = 1;
		//printf("hier%i\n",i);i++;
        ITensor U,S,V,temp;    
        MPS_bins.position(timestep);
        //Print(MPS_bins.A(timestep));
        //Print(MPS_bins.A(timestep+1));
        temp = noprime(MPS_bins.A(timestep)*MPS_bins.A(timestep+1));
        //Print(temp);
        //printf("hier%i\n",i);i++;
        MPS_chain.setA(N,temp);
        //Print(MPS_chain);
        
        //printf("hier%i\n",i);i++;
        //apply the MPO
        MPS *res = new MPS;
		*res = MPS(N);
        MPS_chain = Feedback::apply_nonherm_MPO(MPO_chain, MPS_chain, *res, args);
        delete res;
        //printf("hier%i\n",i);i++;
        if(timestep % 50 == 0) printf("time %.3f of %.3f -- norm_first_site=%.10f -- stagmag=%.10f -- occ_last_site=%.10f -- detector=%.10f\n",
										timestep*dt, t_end, Feedback::norm(MPS_chain.A(1), Site), Feedback::stag_magn(MPS_chain), Feedback::occupation(MPS_chain, N),trans_det);
        if(timestep % 1 == 0)
        {
            fprintf(datastream, "%.3f ", (timestep)*dt);
            for (int i=N;i>0;i--) fprintf(datastream, "\t %.5e ", Feedback::occupation(MPS_chain, i));
            fprintf(datastream, "\t %.10f \n", Feedback::norm(MPS_chain.A(1), Site));
            fprintf(stagmagstream, "%.3f \t %.5e \n", (timestep)*dt, Feedback::stag_magn(MPS_chain));
            fprintf(detectorstream, "%.3f \t %.5e \t %.5e \n", (timestep)*dt, trans_det, dissipated_signal);
        }
        //decompose and write the pair back into the bin-chain (different positions)
        temp = MPS_chain.A(N);   
        //Print(MPS_chain);
        //printf("hier%i\n",i);i++;   
        //U = ITensor(bin[timestep],commonIndex(MPS_bins.A(timestep),MPS_bins.A(timestep-1),Link));
        V = ITensor(s[N],commonIndex(MPS_chain.A(N-1),MPS_chain.A(N),Link),commonIndex(MPS_bins.A(timestep+1),MPS_bins.A(timestep+2)));
        svd(temp,U,S,V);
        //Print(U);
        //Print(V*S);
        //printf("hier%i\n",i);i++;
        MPS_bins.setA(timestep,U);
        MPS_bins.setA(timestep+1,V*S);
        dissipated_signal = Feedback::transient_detector(U*S);
        trans_det += Feedback::transient_detector(U*S);
        //print(MPS_bins);
        //printf("hier%i\n",i);i++;
    }
    
    void Feedback::swap_feedbackbin_to_system(int timestep)
    {
        ITensor U,S,V,temp;
        Index feedbacktimebin = findtype(MPS_bins.A(timestep+1-Nfb),Site);
        for(int i=timestep+1-Nfb; i<timestep; i++)
        {
            temp = MPS_bins.A(i)*MPS_bins.A(i+1);
            //Print(temp);
            U = ITensor(feedbacktimebin, commonIndex(MPS_bins.A(i+1),MPS_bins.A(i+2)));
            svd(temp,U,S,V,args);
            MPS_bins.setA(i,V);
            MPS_bins.setA(i+1,U*S);
        }
    }
    
    void Feedback::swap_feedbackbin_to_origin(int timestep)
    {
        ITensor temp,U,S,V;
        Index feedbacktimebin = findtype(MPS_bins.A(timestep),Site);
        // loop ends when feedback bin is located one bin right of its original position
        // because in the last step, the orthocenter needs to be placed
        // at the future feedback bin for the next time step
        for(int i=timestep; i>timestep-Nfb+2; i--) 
        {
            temp = MPS_bins.A(i-1)*MPS_bins.A(i);
            //Print(temp);
            U = ITensor(feedbacktimebin, leftLinkInd(MPS_bins,i-1));
            svd(temp,U,S,V,args);
            MPS_bins.setA(i,V);
            MPS_bins.setA(i-1,U*S);
        }
        temp = MPS_bins.A(timestep-Nfb+1)*MPS_bins.A(timestep-Nfb+2); 
        //Print(temp);
        U = ITensor(feedbacktimebin, leftLinkInd(MPS_bins,timestep-Nfb+1));
        svd(temp,U,S,V,args);
        MPS_bins.setA(timestep-Nfb+2,V*S);
        MPS_bins.setA(timestep-Nfb+1,U);
        //print("backswap done:\n");
        //Print(MPS_bins);
    }

}//eof namespace




 
