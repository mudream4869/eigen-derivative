#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <Eigen/Dense>

#include "Derivative.cpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using Eigen::Wrapper;
using Eigen::VariableDerivative;

using std::function;
using std::vector;

VectorXd GaussNewtonMethod(vector<Wrapper> fs, VectorXd x){
    int x_size = x.size(), f_size = fs.size();

    vector< vector<Wrapper> > Jac( f_size, vector<Wrapper>(x_size) );
    
    for(int lf = 0;lf < f_size;lf++)
        for(int lx = 0;lx < x_size;lx++)
            Jac[lf][lx] = fs[lf].diffPartial(lx);
    
    for(int iter = 0;iter < 200000;iter++){
        MatrixXd J(f_size, x_size);
        for(int lf = 0;lf < f_size;lf++)
            for(int lx = 0;lx < x_size;lx++)
                J(lf, lx) = Jac[lf][lx](x);

        VectorXd rx(f_size);
        for(int lf = 0;lf < f_size;lf++)
            rx[lf] = fs[lf](x);

        MatrixXd invJ = J.completeOrthogonalDecomposition().pseudoInverse();
        x += invJ*rx;

        std::cout << x.transpose() << std::endl;
    }

    return x;
}

int main(){
    // Solve : xx + yy + zz = 2
    //          x + y + z = 1

    Wrapper x = new VariableDerivative(0),
            y = new VariableDerivative(1),
            z = new VariableDerivative(2);

    Wrapper f1 = x*x + y*y + z*z - 2, f2 = x + y + z - 1;

    VectorXd v(3);
    v << 1, 1, 1;

    std::cout << GaussNewtonMethod({f1, f2}, v).transpose() << std::endl;

    return 0;
}
