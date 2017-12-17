# eigen-derivative
Calculate multivariable function's derivative base on Eigen &lt;3

# Example

```c++
#include <iostream>
#include "Derivative.cpp"

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

    // Divide

    Derivative h = 1/(x*x + 1);
    std::cout << h(v3) << std::endl;
    std::cout << h.diffPartial(0)(v3) << std::endl;

    return 0;
}
```

# TODO 

- [x] ostream output
- [ ] GC
- [ ] more operators
