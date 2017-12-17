#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <Eigen/Dense>

#include "Derivative.cpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using Eigen::Derivative;

using std::function;
using std::vector;

VectorXd LevenbergMarquardt(vector<Derivative> fs, VectorXd x){
    int x_size = x.size(), f_size = fs.size();

    vector< vector<Derivative> > Jac( f_size, vector<Derivative>(x_size) );
    
    for(int lf = 0;lf < f_size;lf++)
        for(int lx = 0;lx < x_size;lx++)
            Jac[lf][lx] = fs[lf].diffPartial(lx);
    
    // Set eps to be very larg
    double mu = 0.01, eps = 10000000;

    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    for(int iter = 0;iter < 200;iter++){
        MatrixXd J(f_size, x_size);
        for(int lf = 0;lf < f_size;lf++)
            for(int lx = 0;lx < x_size;lx++)
                J(lf, lx) = Jac[lf][lx](x);

        VectorXd rx(f_size);
        for(int lf = 0;lf < f_size;lf++)
            rx[lf] = fs[lf](x);

        double rx_norm = rx.norm();

        VectorXd delta = (mu*I + J.transpose()*J).inverse()*J.transpose()*rx;
        x -= delta;

        if(eps < rx_norm)
            mu *= 2;
        else if(eps > rx_norm)
            mu /= 2;

        if(rx_norm < 0.0001 or delta.norm() < 0.0001)
            break;

        eps = rx_norm;
    }

    return x;
}

int main(){
    // Solve : xx + yy + zz = 2
    //          x + y + 2*z = 1

    Derivative x = Derivative::Variable(0),
               y = Derivative::Variable(1),
               z = Derivative::Variable(2);

    Derivative f1 = x*x + y*y + z*z - 2, f2 = x + y + 2*z - 1;

    VectorXd v(3);
    v << 1, 1, 1;

    std::cout << LevenbergMarquardt({f1, f2}, v).transpose() << std::endl;

    return 0;
}
