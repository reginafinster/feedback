#ifndef __feedback__
#define __feedback__

#include "includes.h"

namespace itensor{
class Feedback
    {
    protected:
        double cutoff;
        int maxm;
        int steadybins;
        double dt;
        double t_end;
        int N;
        double J;
        double Jz;
        double Omega;
        double Gamma;
        double phi;
        double inphi;
        double phi_start;
        double phi_stepsize;
        int phi_steps;
        int Nbin;
        int Nfb;
        
        double ttotal;
        
        double INTfak;
        double FBfak;
        double Magn;
		double trans_stagmag;
        double trans_det;
        double dissipated_signal;

        
        int enable_feedback;
        int read_in_shifts;
        int log_transient_stagmag;
        int log_transient_detector;
        
        int initial;
        int init_pos;
        int init_pos1, init_pos2, init_pos3, init_pos4, init_pos5;
        int one_excitation;
        int neel_state_first_up;        
        int neel_state_first_down;        
        int allup_state;   
        int various_excitations;
        
        Args args;
        std::string file_identifier;
        std::string file_identifier1;
        std::string file_identifier2;
        std::string file_identifier3;
        FILE* datastream;
        FILE* stagmagstream;
        FILE* steadystream;
        FILE* detectorstream;
        std::vector<Index> bin;
        std::vector<Index> s;
        MPS MPS_bins;
        MPS MPS_chain;
        ITensor Hsys;
        ITensor Hint;
        ITensor Hfb;
        ITensor Uevo;
        MPO MPO_chain;
        
        ITensor identity(Index& s);
        ITensor MySz(Index& s);
        ITensor MySp(Index& s);
        ITensor MySm(Index& s);
        
        double TLS_occupation(ITensor psitensor, IndexType type);
        double occupation(MPS& MPS_chain, int i);
        double integrate_detector(MPS& MPS_bins,int startbin,int lastbin);
        double transient_detector(ITensor MPSTensor);
        double stag_magn(MPS& MPS_chain);
        double norm(ITensor psi, IndexType type);
        std::string get_date_as_string();
        void setup_dir_and_filename();
        void setup_indices();
        void setup_shifts(double* h);
        void initialize_MPS();
        void set_up_MPO();
        void apply_Uevo_without_feedback(int timestep);
        void apply_Uevo_with_feedback(int timestep);
        void swap_feedbackbin_to_system(int timestep);
        void swap_feedbackbin_to_origin(int timestep);
        MPS apply_nonherm_MPO(MPO& Uevo, MPS& psi, MPS& res, Args args);
        MPO multiply_two_nonherm_MPO(const MPO& MPO_A, const MPO& MPO_B);
    public:
        Feedback();
        Feedback(std::string infile);
        void timeevolve_system();
    };
}
#endif
 
