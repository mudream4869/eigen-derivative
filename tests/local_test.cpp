#include <iostream>
#include "Derivative.h"

using Eigen::Derivative;

int main(){
    Derivative x = Derivative::Variable(0), 
               y = Derivative::Variable(1);
    
    {
        auto f = x*x*x + 3*x - 5*y;
        std::cout << f.diffPartial(0) << std::endl;
    }

    {
        auto f = x*x*y + x*y - 5*y*y;
        std::cout << f << std::endl;
        std::cout << f.diffPartial(1) << std::endl;
    }

    return 0;
}
