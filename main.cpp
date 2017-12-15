#include <iostream>
#include "Derivative.cpp"

int main(){
    // f(x, y) = 3x + 4y
    // dfdx = 3
    // dfdy = 4

    VectorXd v(2);
    v << 3, 4;
    Derivative* dv = new LinearDerivative(v);
    std::cout << dv->call(v) << std::endl;

    
    return 0;
}
