#include <Eigen/Dense>
#include <functional>
#include <vector>

using Eigen::VectorXd;
using std::function;
using std::vector;

typedef function<double(VectorXd)> MultiVar;


class Derivative{
public:
    Derivative(){}
    virtual Derivative diffPartial(int index){}
    virtual double operator()(VectorXd vec){}
};


class SingleDerivative : public Derivative{
private:
    MultiVar inst;
    vector<Derivative*> inst_pd;

public:
    SingleDerivative(MultiVar _inst, vector<Derivative*> _inst_pd):
        inst(_inst),inst_pd(_inst_pd){}

    Derivative diffPartial(int index){
        return *inst_pd[index];
    }

    double operator()(VectorXd vec){
        return inst(vec);
    }
};


class ZeroDerivative : public Derivative{
public:
    ZeroDerivative(){}
    
    Derivative diffPartial(int index){
        return *this;
    }
    
    double operator()(VectorXd vec){
        return 0;
    }
};



class DerivativeAdd : public Derivative{
private:
    Derivative a, b;

public:
    DerivativeAdd(Derivative _a, Derivative _b){a = _a, b = _b;}
    Derivative diffPartial(int index);
    double operator()(VectorXd vec){
        return a(vec) + b(vec);
    }
};


Derivative operator+(Derivative a, Derivative b){return DerivativeAdd(a, b);}
