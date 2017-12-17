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

typedef function<double(VectorXd)> FuncDV;
typedef function<VectorXd(VectorXd)> FuncVV;
typedef function<MatrixXd(VectorXd)> FuncMV;

MatrixXd MakePositiveSemidefinite(MatrixXd mat, int dim){
    double left = 0, right = std::max(-mat.minCoeff(), .0) + 1;
    while(right - left >= 1e-4){
        double mid = (left+right)/2;
        auto test = mat + mid*MatrixXd::Identity(dim, dim);
        // From https://stackoverflow.com/questions/35227131/eigen-check-if-matrix-is-positive-semi-definite

        Eigen::LLT<Eigen::MatrixXd> lltOfA(test);
        if(lltOfA.info() == Eigen::NumericalIssue)
            left = mid;
        else
            right = mid;
    }
    return mat + right*MatrixXd::Identity(dim, dim);
}

VectorXd doIPM(
    FuncDV F, FuncVV DelF, FuncMV LaplaceF,
    FuncVV H, FuncMV DelH, vector<FuncMV> LaplaceH,
    int dimX, int dimH, VectorXd initX
){
    double mu = 0.0001;

    std::function<bool(VectorXd)> feasible = [dimH](VectorXd _f){
        for(int lx = 0;lx < dimH;lx++)
            if(_f[lx] < 0)
                return false;
        return true; 
    };

    std::function<double(VectorXd, VectorXd, VectorXd)> 
        L = [mu, F, H, dimH](VectorXd x, VectorXd y, VectorXd w){
        double barrier = 0;
        for(int lx = 0;lx < dimH;lx++)
            barrier += log(w[lx]);
        return F(x) - mu*barrier - y.transpose()*(H(x) - w);
    };

    VectorXd x = initX, y, w = H(initX);

    { 
        MatrixXd A = DelH(x).transpose();
        MatrixXd invA = A.completeOrthogonalDecomposition().pseudoInverse();
        y = invA*DelF(x);
    }
    
    for(;;){
        VectorXd h = H(x), e = VectorXd::Ones(dimH);
        MatrixXd W = w.asDiagonal(), invY = y.asDiagonal().inverse();

        MatrixXd Hess = LaplaceF(x);
        for(int lx = 0;lx < dimH;lx++)
            Hess -= y[lx]*LaplaceH[lx](x);

        // Find ~H = H + lambda * I
        Hess = MakePositiveSemidefinite(Hess, dimX);

        MatrixXd A = DelH(x);
        MatrixXd M1(dimX + dimH, dimX + dimH);

        M1 << -Hess, A.transpose(),
                  A, W*invY       ;

        VectorXd V1(dimX + dimH);
        
        V1 << DelF(x) - A.transpose()*y,
                  - h + mu*invY*e      ;

        VectorXd dxy = M1.inverse()*V1;
        VectorXd dx, dy, dw;

        dx = dxy.head(dimX);
        dy = dxy.tail(dimH);
        dw = invY*mu*e - W*e - invY*W*dy;

        x += dx, y += dy, w += dw;

        if(dx.norm() <= 0.001)
            break;
    }

    return x;
}

// Solve : min  obj_f
//         sub  con_hs >= 0 

VectorXd Derivative_IPM(Derivative obj_f, vector<Derivative> con_hs, VectorXd start_guess){
    int x_size = start_guess.size(), h_size = con_hs.size();

    vector<Derivative> v1w(x_size);
    vector< vector<Derivative> > v2w(x_size, v1w);

    vector<Derivative> f_gradient = v1w;
    vector< vector<Derivative> > f_hess = v2w;

    vector< vector<Derivative> > hs_gradient(h_size, v1w);
    vector< vector< vector<Derivative> > > hs_hess(h_size, v2w);

    // Prebuild differiential function
    
    for(int lx = 0;lx < x_size;lx++)
        f_gradient[lx] = obj_f.diffPartial(lx);

    for(int lx = 0;lx < x_size;lx++)
        for(int ly = 0;ly < x_size;ly++)
            f_hess[lx][ly] = f_gradient[lx].diffPartial(ly);
    
    for(int lh = 0;lh < h_size;lh++){
        for(int lx = 0;lx < x_size;lx++)
            hs_gradient[lh][lx] = con_hs[lh].diffPartial(lx);

        for(int lx = 0;lx < x_size;lx++)
            for(int ly = 0;ly < x_size;ly++)
                hs_hess[lh][lx][ly] = hs_gradient[lh][lx].diffPartial(ly);
    }

    FuncDV F = [obj_f](VectorXd x){
        return obj_f(x);
    };

    FuncVV DelF = [f_gradient, x_size](VectorXd x){
        VectorXd ret(x_size);
        for(int lx = 0;lx < x_size;lx++)
            ret[lx] = f_gradient[lx](x);
        return ret;
    };

    FuncMV LaplaceF = [f_hess, x_size](VectorXd x){ 
        MatrixXd ret(x_size, x_size);
        for(int lx = 0;lx < x_size;lx++)
            for(int ly = 0;ly < x_size;ly++)
                ret(lx, ly) = f_hess[lx][ly](x);
        return ret; 
    };

    FuncVV H = [con_hs, h_size](VectorXd x){
        VectorXd ret(h_size);
        for(int lx = 0;lx < h_size;lx++)
            ret[lx] = con_hs[lx](x);
        return ret;
    };

    FuncMV DelH = [hs_gradient, h_size, x_size](VectorXd x){
        MatrixXd A(h_size, x_size);
        for(int lh = 0;lh < h_size;lh++)
            for(int lx = 0;lx < x_size;lx++)
                A(lh, lx) = hs_gradient[lh][lx](x);
        return A;
    };
    
    vector<FuncMV> LaplaceH(h_size);
    
    for(int lh = 0;lh < h_size;lh++){
        LaplaceH[lh] = [hs_hess, lh, x_size](VectorXd x){
            MatrixXd ret(x_size, x_size);
            for(int lx = 0;lx < x_size;lx++)
                for(int ly = 0;ly < x_size;ly++)
                    ret(lx, ly) = hs_hess[lh][lx][ly](x);
            return ret;
        };
    };
 
    return doIPM(F, DelF, LaplaceF, H, DelH, LaplaceH, start_guess.size(), con_hs.size(), start_guess); 
}

int main(){
    // Reference : http://www.princeton.edu/~rvdb/tex/talks/MLSS_LaPalma/LaPalma3.pdf
    // Solve : min  x+y
    //         sub  xx + yy >= 1 
    //              x >= 0 
    //              y >= 0 

    Derivative par_x = Derivative::Variable(0), par_y = Derivative::Variable(1);
    Derivative obj_f = par_x + par_y;

    Derivative con_h1 = par_x*par_x + par_y*par_y - 1,
               con_h2 = par_x,
               con_h3 = par_y;

    // Multiple initial value
    VectorXd x(2);
    
    x << 1, 1;
    std::cout << "x initial as " << x.transpose() << std::endl;
    std::cout << Derivative_IPM(obj_f, {con_h1, con_h2, con_h3}, x).transpose() << std::endl;

    x << 1, 2;
    std::cout << "x initial as " << x.transpose() << std::endl;
    std::cout << Derivative_IPM(obj_f, {con_h1, con_h2, con_h3}, x).transpose() << std::endl;

    x << 2, 1;
    std::cout << "x initial as " << x.transpose() << std::endl;
    std::cout << Derivative_IPM(obj_f, {con_h1, con_h2, con_h3}, x).transpose() << std::endl;
    
    x << 4, -1;
    std::cout << "x initial as " << x.transpose() << std::endl;
    std::cout << Derivative_IPM(obj_f, {con_h1, con_h2, con_h3}, x).transpose() << std::endl;

    x << -1, 4;
    std::cout << "x initial as " << x.transpose() << std::endl;
    std::cout << Derivative_IPM(obj_f, {con_h1, con_h2, con_h3}, x).transpose() << std::endl;

    x << 0.2, 0.7;
    std::cout << "x initial as " << x.transpose() << std::endl;
    std::cout << Derivative_IPM(obj_f, {con_h1, con_h2, con_h3}, x).transpose() << std::endl;

    x << 0.5, 0.5;
    std::cout << "x initial as " << x.transpose() << std::endl;
    std::cout << Derivative_IPM(obj_f, {con_h1, con_h2, con_h3}, x).transpose() << std::endl;
    return 0;
}
