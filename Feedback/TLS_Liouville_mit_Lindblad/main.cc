#include "TLS_Liouville.h"

using namespace itensor;

int main(int argc, char* argv[])
{
    if(argc != 2) 
    { 
    printf("Please specify parameter input file."); 
    return 0; 
    }
    
    TLS_Liouville TLS_Liouville(argv[1]);
    TLS_Liouville.timeevolve_system();
    return 0;
}
