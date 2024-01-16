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
        double GammaOut;
        double GammaIn;
        double Kappa;
        double phi;
        int Nfb;
        double ttotal;
        int save_nth_timestep;
        int print_nth_timestep;
        
        Args args;
        std::string file_identifier_timeevo;
        FILE* datastream_timeevo;
        
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
        double excited_occupation(ITensor psi);
        double full_exc_occ(ITensor psi);
        double full_gr_occ(ITensor psi);
        double ground_occupation(ITensor psi);
        double rhosquare(ITensor psi);
        double single_tensor_trace(ITensor psitensor);
        double full_trace();
        ITensor contract_psi(MPS psi, int ttotal);
        double trace_contracted(ITensor psi_contracted);
        double tr_ex_occ(ITensor psi_contracted);
        double tr_gr_occ(ITensor psi_contracted);
        ITensor apply_MPO_for_readout(MPO MPO_Operator);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        ITensor identity(Index s);
        void setup_indexbin();
        void initialize_MPS_with_feedback();
        void initialize_MPS_without_feedback();
        void set_up_Uevo_with_feedback();
        void set_up_Uevo_without_feedback();
        void apply_Uevo_with_feedback(int timestep);
        void apply_Uevo_without_feedback(int timestep);
        void swap_feedbackbin_to_origin(int timestep);
        void swap_feedbackbin_to_system(int timestep);

    public:
        TLS_Liouville();
        TLS_Liouville(std::string infile);
        void timeevolve_system_with_feedback();
        void timeevolve_system_without_feedback();
    };
}
#endif
