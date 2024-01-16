#ifndef __hom__
#define __hom__

#include "includes.h"

namespace itensor{
class HOM
    {
    protected:
        double cutoff;
        int maxm;
        double dt;
        double t_end;
        double Gamma1;
        double Gamma2;
        double Gamma2Max;
        int GammaSteps;
        double phi1;
        double phi2;
        int Nbin;
        int Nfb1;
        int Nfb2;
        double ttotal1;
        double ttotal2;
        int save_every_nth_step;
        int feedback_on_1;
        int feedback_on_2;
        std::string output_dir;
        
        Args args;
        std::string file_identifier1;
        std::string file_identifier2;
        std::string file_identifier3;
        std::string file_identifier4;
        std::string file_identifier6;
        FILE* datastream1;
        FILE* datastream2;
        FILE* datastream3;
        FILE* datastream4;
        FILE* datastream6;
        std::vector<Index> bin1;
        std::vector<Index> bin2;
        MPS psi1;
        MPS psi2;
        ITensor Hint1;
        ITensor Hfb1;
        ITensor Hint2;
        ITensor Hfb2;
        ITensor Uevo1;
        ITensor Uevo2;
        
        double TLS_occupation(ITensor psi, IndexType type);
        double norm(ITensor psi, IndexType type);
        std::string get_date_as_string();
        ITensor identity(Index s);
        void setup_indexbin1();
        void setup_indexbin2();
        void initialize_MPS1_without_feedback();
        void initialize_MPS2_without_feedback();
        void set_up_MPO1_without_feedback();
        void set_up_MPO2_without_feedback();
        void apply_Uevo1_without_feedback(int timestep);
        void apply_Uevo2_without_feedback(int timestep);
        
        double temp_bath_occupation1;
        double temp_bath_occupation2;
        double norm1;
        double norm2;
        
        double g2_ref_sum;
        double g2_ref_sum_norm;
        ITensor annihilator(Index s);
        ITensor creator(Index s);
        double calc_bath_pop(ITensor temp, ITensor Sm_RefBin);
        std::vector<std::vector<double>> g2;
        std::vector<double> pop1;
        std::vector<double> pop2;
        ITensor Sm_RefBin1;
        ITensor Sm_RefBin2;
        ITensor Sp_RefBin1;
        ITensor Sp_RefBin2;
        
        void initialize_MPS1_with_feedback();
        void initialize_MPS2_with_feedback();
        void set_up_MPO1_with_feedback();
        void set_up_MPO2_with_feedback();
        void apply_Uevo1_with_feedback(int timestep);
        void apply_Uevo2_with_feedback(int timestep);
        void swap_feedbackbin_to_system1(int timestep);
        void swap_feedbackbin_to_system2(int timestep);
        void swap_feedbackbin_to_origin1(int timestep);
        void swap_feedbackbin_to_origin2(int timestep);
    public:
        HOM();
        HOM(std::string infile);
        void timeevolve_system1_with_feedback();
        void timeevolve_system2_with_feedback();
        void timeevolve_system1_without_feedback();
        void timeevolve_system2_without_feedback();
        void calc_g2_without_feedback();
        void calc_g2_with_feedback();
        void loop_over_Gamma();
        void setup_dir_and_filename();
    };
}
#endif
