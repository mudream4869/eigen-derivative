#include <iostream>
#include <chrono>
#include "Derivative.cpp"

using Eigen::Derivative;

struct Timer{
    std::chrono::steady_clock::time_point t1;
    Timer(){
        t1 = std::chrono::steady_clock::now();
    }

    double operator()(){
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    }
};
    

int main(){

    // Variable Derivative

    const int N = 100;
    
    Derivative p = 1;
    for(int lx = 0;lx < N;lx++)
        p = p*Derivative::Variable(lx);

    VectorXd v = VectorXd::Random(N);

    Timer t1; 
    auto dp = p.diffPartial(0);
    std::cout << "Partial differential cost " << t1() << "s" << std::endl;
    std::cout << dp << std::endl;

    Timer t2;
    std::cout << dp(v) << std::endl; 
    std::cout << "Calculate value cost " << t2() << "s" << std::endl;

    Timer t3;
    for(int lx = 0;lx < N;lx++)
        p = p.diffPartial(lx);
    std::cout << "Partial diff all variables cost " << t3() << "s" << std::endl;

    return 0;
}
