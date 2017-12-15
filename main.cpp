#include <iostream>
#include "Derivative.cpp"

int main(){
    // f(x, y) = 3x + 4y
    // dfdx = 3
    // dfdy = 4

    VectorXd v(2);
    v << 3, 4;
    Wrapper dv = new LinearDerivative(v);
    std::cout << (dv + dv).diffPartial(0)(v) << std::endl;;

    // f(x, y, z) = 3x + 4y + 5zz
    // dfdx = 3, dfdy = 4, dfdz = 10z

    VectorXd v1(3), v2(3);
    v1 << 3, 4, 0;
    v2 << 0, 0, 1;
    
    Wrapper five = new ConstantDerivative(5);
    Wrapper xpy = new LinearDerivative(v1);
    Wrapper z = new LinearDerivative(v2);
    Wrapper f = xpy + five*z*z;

    std::cout << f.diffPartial(0)(v1) << std::endl;
    std::cout << f.diffPartial(2)(v2) << std::endl;

    // Variable Derivative

    VectorXd v3(2);
    v3 << 3, 4;
    Wrapper x = new VariableDerivative(0);
    Wrapper y = new VariableDerivative(1);
    Wrapper g = x*x + x*y + y*y;
    std::cout << g.diffPartial(0)(v) << std::endl;
    std::cout << g.diffPartial(1)(v) << std::endl;

    // Divide

    Wrapper h = 1/(x*x + 1);
    std::cout << h(v3) << std::endl;
    std::cout << h.diffPartial(0)(v3) << std::endl;

    return 0;
}
