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
        int print_every_nth_step;
        int save_every_nth_step;
        double ttotal;
        double init_22;
        double init_11;
        double benchmark_factor;
        int spectrum_steps;
        double spectrum_interval;
        
        Args args;
        std::string file_identifier_timeevo;
        std::string file_identifier_correlation;
        std::string file_identifier_spectrum;
        FILE* datastream_timeevo;
        FILE* datastream_correlation;
        FILE* datastream_spectrum;
        std::vector<Index> bin;
        MPS psi;
        ITensor Hsys;
        ITensor Hint;
        ITensor Hfb;
        ITensor Uevo;
        
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
        
        double occupation(ITensor psi, IndexType type);
        double ground_occupation(ITensor psi, IndexType type);
        double norm(ITensor psi);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        ITensor identity(Index s);
        void setup_indexbin();
        void calc_exact_variables(int timestep);
        void initialize_MPS();
        void set_up_MPO();
        void apply_Uevo(int timestep);
        void swap_feedbackbin_to_system(int timestep);
        void swap_feedbackbin_to_origin(int timestep);
        ITensor annihilator(Index index);
        std::complex<double> calc_g1_bath(Index timebin, ITensor temp);
        std::complex<double> calc_g2_bath(Index timebin, ITensor temp);
        std::complex<double> calc_pop_bath(ITensor temp);
        
    public:
        Mollow();
        Mollow(std::string infile);
        double get_ttotal() const { return ttotal; }
        double get_spectrum_steps() const { return spectrum_steps; }
        void timeevolve_system();
        void calc_correlation();
        void calc_spectrum();
    };
}
#endif
