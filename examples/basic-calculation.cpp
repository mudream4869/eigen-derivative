#include <iostream>
#include "Derivative.cpp"

using Eigen::Wrapper;

int main(){

    // Variable Derivative

    VectorXd v3(2);
    v3 << 3, 4;
    Wrapper x = Wrapper::Variable(0);
    Wrapper y = Wrapper::Variable(1);
    Wrapper g = x*x + x*y + y*y;
    std::cout << g.diffPartial(0)(v3) << std::endl;
    std::cout << g.diffPartial(1)(v3) << std::endl;

    // Divide

    Wrapper h = 1/(x*x + 1);
    std::cout << h(v3) << std::endl;
    std::cout << h.diffPartial(0)(v3) << std::endl;

    return 0;
}
