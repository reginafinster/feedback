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
        double LindbladGamma;
        double phi;
        int Nfb;
        double ttotal;
        Args args;
        std::string file_identifier;
        std::string init_cond;
        std::string use_feedback;
        FILE* datastream;
        std::vector<Index> bin;
        MPS psi;
        ITensor Lsys;
        ITensor Lint;
        ITensor Lfb;
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
        void initialize_MPS_without_feedback();
        void set_up_Uevo_without_feedback();
        void apply_Uevo_without_feedback(int timestep);
        void initialize_MPS_with_feedback();
        void set_up_Uevo_with_feedback();
        void apply_Uevo_with_feedback(int timestep);
        void swap_feedbackbin_to_origin(int timestep);
        void swap_feedbackbin_to_system(int timestep);

    public:
        TLS_Liouville();
        TLS_Liouville(std::string infile);
        void timeevolve_system_without_feedback();
        void timeevolve_system_with_feedback();
    };
}
#endif
