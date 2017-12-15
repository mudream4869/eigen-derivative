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
    virtual Derivative* diffPartial(int index){}
    virtual double call(VectorXd vec){}
};


class SingleDerivative : public Derivative{
private:
    MultiVar inst;
    vector<Derivative*> inst_pd;

public:
    SingleDerivative(MultiVar _inst, vector<Derivative*> _inst_pd):
        inst(_inst),inst_pd(_inst_pd){}

    Derivative* diffPartial(int index){
        return inst_pd[index];
    }

    double call(VectorXd vec){
        return inst(vec);
    }
};


class ConstantDerivative : public Derivative{
private:
    double a;

public:
    ConstantDerivative(double _a):a(_a){}
    
    Derivative* diffPartial(int index){
        return new ConstantDerivative(0);
    }
    
    double call(VectorXd vec){
        return a;
    }
};


class LinearDerivative : public Derivative{
private:
    VectorXd v;

public:
    LinearDerivative(VectorXd _v){
        v = _v;
    }
 
    Derivative* diffPartial(int index){
        return new ConstantDerivative(v[index]);
    }
    
    double call(VectorXd vec){
        return v.transpose()*vec;
    }
};


class DerivativeAdd : public Derivative{
private:
    Derivative* a;
    Derivative* b;

public:
    DerivativeAdd(Derivative* _a, Derivative* _b){a = _a, b = _b;}
    Derivative* diffPartial(int index);
    double operator()(VectorXd vec){
        return a->call(vec) + b->call(vec);
    }
};


class DerivativeMultiply : public Derivative{
private:
    Derivative* a;
    Derivative* b;

public:
    DerivativeMultiply(Derivative* _a, Derivative* _b){a = _a, b = _b;}
    Derivative* diffPartial(int index);
    double operator()(VectorXd vec){
        return a->call(vec) * b->call(vec);
    }
};


class Wrapper{
public:
    Derivative* inst;
    Wrapper(Derivative* _inst):inst(_inst){}
    Wrapper diffPartial(int index){return inst->diffPartial(index);}
    double operator()(VectorXd vec){return inst->call(vec);}
};


Wrapper operator+(Wrapper a, Wrapper b){
    return new DerivativeAdd(a.inst, b.inst);
}


Wrapper operator*(Wrapper a, Wrapper b){
    return new DerivativeMultiply(a.inst, b.inst);
}
