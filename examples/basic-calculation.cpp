#include <iostream>
#include "Derivative.h"

using Eigen::VectorXd;
using Eigen::Derivative;

int main(){

    // Variable Derivative

    VectorXd v3(2);
    v3 << 3, 4;
    Derivative x = Derivative::Variable(0);
    Derivative y = Derivative::Variable(1);
    Derivative g = x*x + x*y + y*y;
    std::cout << g.diffPartial(0)(v3) << std::endl;
    std::cout << g.diffPartial(1)(v3) << std::endl;

    std::cout << g.diffPartial(0) << std::endl;

    // Divide

    Derivative h = 1/(x*x + 1);
    std::cout << h(v3) << std::endl;
    std::cout << h.diffPartial(0)(v3) << std::endl;

    // More
    
    Derivative p = exp(x*x) + log(x);
    std::cout << p << std::endl;
    std::cout << p.diffPartial(0) << std::endl;

    Derivative p2 = pow(x, 2.5);
    std::cout << p2 << std::endl;
    std::cout << p2.diffPartial(0).diffPartial(0) << std::endl;

    return 0;
}
