#include <iostream>
#include "Derivative.cpp"

using Eigen::Derivative;
using Eigen::Matrix;

typedef Matrix<Derivative, 2, 1> M2F;
typedef Matrix<Derivative, 2, 2> M2x2F;

int main(){
    auto x = Derivative::Variable(0),
         y = Derivative::Variable(1);
    auto f = x + y + x*y + log(x);

    M2F df1, df2, df3;
    df1 << f.diffPartial(0), f.diffPartial(1);
    df2 << f.diffPartial(1), f.diffPartial(0);
    
    df3 = df1 + df2;
    M2x2F m2x2f = df1*df2.transpose();
    Derivative a = df1.dot(df2), b = df1.sum();

    std::cout << a << std::endl;
    std::cout << b << std::endl;

    return 0;
}
