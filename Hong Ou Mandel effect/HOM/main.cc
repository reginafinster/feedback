#include "hom.h"

using namespace itensor;

int main(int argc, char* argv[])
{
    if(argc != 2) 
    { 
    printf("Please specify parameter input file."); 
    return 0; 
    }
    
    HOM hom(argv[1]);
    //hom.setup_dir_and_filename();
    //hom.timeevolve_system1_with_feedback();
    //hom.timeevolve_system2_with_feedback();
    //hom.calc_g2_with_feedback_on_both_systems();
    hom.loop_over_Gamma();
    return 0;
}
