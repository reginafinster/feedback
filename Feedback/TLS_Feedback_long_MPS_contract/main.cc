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
    TLS_feedback.timeevolve_system_with_feedback();
    std::cout << "main" << std::endl;
    return 0;
    
}
