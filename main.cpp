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

    std::cout << f(v1) << std::endl;


    
    return 0;
}
