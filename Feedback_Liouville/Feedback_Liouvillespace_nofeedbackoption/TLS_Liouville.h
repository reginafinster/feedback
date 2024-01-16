#ifndef __TLS_Liouville__
#define __TLS_Liouville__

#include "includes.h"

namespace itensor{
class TLS_Liouville
    {
    protected:
        double cutoff;
        int maxm;
        double dt;
        double t_end;
        double Omega;
        double Gamma;
        double LindbladLossGamma;
        double LindbladPumpGamma;
        double GammaPureDephasing;
        double phi;
        int Nfb;
        double ttotal;
        int save_nth_timestep;
        int print_nth_timestep;
        double spectrum_steps;
        double spectrum_interval;  
        
        
        int referencetimebin;
        ITensor Sm_RefBin;
        std::vector<Complex> g1_bath;
        std::vector<Complex> g2_bath;
        std::vector<Complex> pop_bath;
        std::vector<Complex> spectrum;
        ITensor g1_tensor;
        ITensor g2_tensor;
        ITensor pop_bath_tensor;
        double dw;
        double om;
        
        Args args;
        std::string file_identifier_timeevo;
        std::string file_identifier_correlation;
        std::string file_identifier_spectrum;
        FILE* datastream_timeevo;
        FILE* datastream_correlation;
        FILE* datastream_spectrum;
        
        std::string init_cond;
        std::string use_feedback;
        FILE* datastream;
        std::vector<Index> bin;
        MPS psi;
        ITensor Lsys;
        ITensor Lint;
        ITensor Lfb;
        ITensor LindbladLoss;
        ITensor LindbladPump;
        ITensor PureDephasing;
        ITensor Lindblad;
        ITensor Uevo;
        double excited_TLS_occupation(ITensor psi, IndexType type);
        double ground_TLS_occupation(ITensor psi, IndexType type);
        double rhosquare(ITensor psi, IndexType type);
        double trace(ITensor psitensor, IndexType type);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        ITensor identity(Index s);
        void setup_indexbin();
        void initialize_MPS_with_feedback();
        void set_up_Uevo_with_feedback();
        void apply_Uevo_with_feedback(int timestep);
        void swap_feedbackbin_to_origin(int timestep);
        void swap_feedbackbin_to_system(int timestep);
        
        ITensor annihilator(Index index);
        std::complex<double> calc_g1_bath(Index timebin, ITensor temp);
        std::complex<double> calc_g2_bath(Index timebin, ITensor temp);
        std::complex<double> calc_pop_bath(ITensor temp);

    public:
        TLS_Liouville();
        TLS_Liouville(std::string infile);
        void timeevolve_system_with_feedback();
        void calc_correlation();
        void calc_spectrum();
    };
}
#endif
