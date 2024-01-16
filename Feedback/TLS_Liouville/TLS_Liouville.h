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
        double phi;
        int Nfb;
        double ttotal;
        
        Args args;
        std::string file_identifier;
        std::string init_cond;
        FILE* datastream;
        std::vector<Index> bin;
        MPS psi;
        ITensor Lsys;
        ITensor Lint;
        ITensor Uevo;
        
        double TLS_occupation(ITensor psi, IndexType type);
        double norm(ITensor psi, IndexType type);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        ITensor identity(Index s);
        void setup_indexbin();
        void initialize_MPS();
        void set_up_operator();
        void apply_Levo(int timestep);
    public:
        TLS_Liouville();
        TLS_Liouville(std::string infile);
        void timeevolve_system();
    };
}
#endif
