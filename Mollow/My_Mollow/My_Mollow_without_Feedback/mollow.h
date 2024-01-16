#ifndef __Mollow__
#define __Mollow__

#include "includes.h"

namespace itensor{
class Mollow
    {
    protected:
        double cutoff;
        int maxm;
        double dt;
        double t_end;
        double Omega;
        double Gamma;
        double phi;
        int Nbin;
        int Nfb;
        double ttotal;
        double init_EE;
        double init_GG;
        double benchmark_factor;
        
        Args args;
        std::string file_identifier;
        std::string init_cond;
        FILE* datastream;
        std::vector<Index> bin;
        MPS psi;
        ITensor Hsys;
        ITensor Hint;
        ITensor Hfb;
        ITensor Uevo;
        
        double tGamma;
        double tOmega;
        double bOm;
        double bGam;    
        double exact;
        double eG;
        double sOm;
        double cOm;
        double S0;
        double W0;
        
        double occupation(ITensor psi, IndexType type);
        double ground_occupation(ITensor psi, IndexType type);
        double coherence(ITensor psitensor, IndexType type);
        double norm(ITensor psi);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        ITensor identity(Index s);
        void setup_indexbin();
        void initialize_MPS_without_feedback();
        void set_up_MPO_without_feedback();
        void apply_Uevo_without_feedback(int timestep);
        void calc_exact_variables(int timestep);
//         void set_up_MPO_with_feedback();
//         void apply_Uevo_with_feedback(int timestep);
//         void swap_feedbackbin_to_system(int timestep);
//         void swap_feedbackbin_to_origin(int timestep);
    public:
        Mollow();
        Mollow(std::string infile);
        void timeevolve_system_with_feedback();
        void timeevolve_system_without_feedback();
    };
}
#endif
