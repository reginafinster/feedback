#include "biexciton.h"

using namespace itensor;

int main(int argc, char* argv[])
{
    if(argc != 2) 
    { 
    printf("Please specify parameter input file."); 
    return 0; 
    }
    
    Biexciton my_biexciton(argv[1]);
    
    if(my_biexciton.get_spectrum_steps() > my_biexciton.get_ttotal())
    {
        println("EXITING: Spectrum steps must be smaller than t_end/dt");
        return 0;
    }
    my_biexciton.timeevolve_system();
    my_biexciton.calc_correlation();
    my_biexciton.calc_spectrum();
    return 0;
}
