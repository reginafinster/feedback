#ifndef __feedback__
#define __feedback__

#include "includes.h"

namespace itensor{
class Feedback
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
        
        double INTfak;
        double FBfak;
        
        
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
        
        double TLS_occupation(ITensor psi, IndexType type);
        double norm(ITensor psi, IndexType type);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        ITensor identity(Index s);
        void setup_indexbin();
        void initialize_MPS();
        void set_up_MPO_without_feedback();
        void apply_Uevo_without_feedback(int timestep);
        void set_up_MPO_with_feedback();
        void apply_Uevo_with_feedback(int timestep);
        void swap_feedbackbin_to_system(int timestep);
        void swap_feedbackbin_to_origin(int timestep);
    public:
        Feedback();
        Feedback(std::string infile);
        void timeevolve_system_with_feedback();
        void timeevolve_system_without_feedback();
    };
}
#endif
