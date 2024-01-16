#include "feedback.h"

using namespace itensor;

int main(int argc, char* argv[])
{
    if(argc != 2) 
    { 
    printf("Please specify parameter input file."); 
    return 0; 
    }
    
    Feedback TLS_feedback(argv[1]);
    // TLS_feedback.timeevolve_system_without_feedback();
    // TLS_feedback.calc_g2_without_feedback();
    TLS_feedback.loop_over_Gamma();
    return 0;
}
