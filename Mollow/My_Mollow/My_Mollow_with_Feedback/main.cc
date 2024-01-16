#include "mollow.h"

using namespace itensor;

int main(int argc, char* argv[])
{
    if(argc != 2) 
    { 
    printf("Please specify parameter input file."); 
    return 0; 
    }
    
    Mollow mollow(argv[1]);
    
    if(mollow.get_spectrum_steps() > mollow.get_ttotal())
    {
        println("EXITING: Spectrum steps must be smaller than t_end/dt");
        return 0;
    }
    
    mollow.timeevolve_system();
    mollow.calc_correlation();
    //mollow.calc_spectrum();
    
    return 0;
}
