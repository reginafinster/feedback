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
    InputGroup input = InputGroup(argv[1],"input");
    std::string use_feedback = input.getString("use_feedback");
    if (use_feedback=="no")
    {
        TLS_Liouville.timeevolve_system_without_feedback();
    }
    else if (use_feedback=="yes")
    {
        TLS_Liouville.timeevolve_system_with_feedback();
    }
    return 0;
}
