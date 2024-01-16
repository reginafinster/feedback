#include "spinchain_feedback.h"

using namespace itensor;

int main(int argc, char* argv[])
{
    if(argc != 2) 
    { 
    printf("Please specify parameter input file."); 
    return 0; 
    }
    Feedback Spinchain_feedback(argv[1]);
    Spinchain_feedback.timeevolve_system();
    std::cout << "done" << std::endl;
    return 0;
}
