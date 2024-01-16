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
    mollow.timeevolve_system_without_feedback();
    return 0;
}
