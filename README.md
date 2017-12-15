# eigen-derivative
Calculate multivariate derivative function base on Eigen &lt;3

# Example

```c++
#include <iostream>
#include "Derivative.cpp"

int main(){
    VectorXd v(2);
    v << 3, 4;
    Wrapper x = new VariableDerivative(0);
    Wrapper y = new VariableDerivative(1);
    Wrapper g = x*x + x*y + y*y;
    std::cout << g(v) << std::endl;
    std::cout << g.diffPartial(1)(v) << std::endl;

    Wrapper h = 1/(x*x + 1);
    std::cout << h(v) << std::endl;
    std::cout << h.diffPartial(0)(v) << std::endl;

    return 0;
}
```
