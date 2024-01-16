#ifndef __Biexction__
#define __Biexction__

#include "includes.h"

namespace itensor{
class Biexciton
    {
    protected:
        double cutoff;
        int maxm;
        double dt;
        double t_end;
        double Omega_H;
        double Delta;
        double Gamma;
        double phi;
        int Nbin;
        int Nfb;
        int print_every_nth_step;
        int save_every_nth_step;
        double ttotal;
        double init_GG, init_VV, init_HH, init_BB;
        double benchmark_factor;
        int spectrum_steps;
        double spectrum_interval;
        
        Args args;
        std::string file_identifier_timeevo, file_identifier_correlationAforH, file_identifier_correlationBforV, file_identifier_spectrumAforH, file_identifier_spectrumBforV;
        FILE* datastream_timeevo;
        FILE* datastream_correlationAforH;
        FILE* datastream_correlationBforV;
        FILE* datastream_spectrumAforH;
        FILE* datastream_spectrumBforV;
        std::vector<Index> binA, binB;
        MPS psi;
        ITensor Hsys;
        ITensor HintA;
        ITensor HintB;
        ITensor Hint;
        ITensor Hfb;
        ITensor Uevo;
        
        int referencetimebin;
        ITensor Sm_RefBinA, Sm_RefBinB;
        std::vector<Complex> g1_bathA, g2_bathA, pop_bathA, g1_bathB, g2_bathB, pop_bathB, spectrumA, spectrumB;
        double dw, om;
        
        double occupation_system(ITensor psi, int system_level);
        double sum_probabilities(ITensor psitensor);
        double norm(ITensor psi);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        ITensor identity(Index s);
        void setup_indexbins();
        void initialize_MPS();
        void set_up_MPO();
        void apply_Uevo(int timestep);
        void swap_feedbackbin_to_system(int timestep);
        void swap_feedbackbin_to_origin(int timestep);
        ITensor annihilator(Index index);
        std::complex<double> calc_g1_bathA(Index timebinA, Index timebinB, ITensor temp);
        std::complex<double> calc_g2_bathA(Index timebinA, Index timebinB, ITensor temp);
        std::complex<double> calc_pop_bathA(ITensor temp);
        std::complex<double> calc_g1_bathB(Index timebinB, Index timebinA, ITensor temp);
        std::complex<double> calc_g2_bathB(Index timebinB, Index timebinA, ITensor temp);
        std::complex<double> calc_pop_bathB(ITensor temp);
        
    public:
        Biexciton();
        Biexciton(std::string infile);
        double get_ttotal() const { return ttotal; }
        double get_spectrum_steps() const { return spectrum_steps; }
        void timeevolve_system();
        void calc_correlation();
        void calc_spectrum();
    };
}
#endif
