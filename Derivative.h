#include <Eigen/Dense>
#include <functional>
#include <vector>

using Eigen::VectorXd;
using std::function;
using std::vector;

typedef function<double(VectorXd)> MultiVar;

namespace Eigen{

class Derivative{
public:
    Derivative(){}
    virtual Derivative* diffPartial(int index){}
    virtual double call(VectorXd vec) const {}
};


class SingleDerivative : public Derivative{
private:
    MultiVar inst;
    function<Derivative*(int)> inst_pd;

public:
    SingleDerivative(MultiVar _inst, function<Derivative*(int)> _inst_pd):
        inst(_inst),inst_pd(_inst_pd){}

    Derivative* diffPartial(int index){
        return inst_pd(index);
    }

    double call(VectorXd vec) const {
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
    
    double call(VectorXd vec) const {
        return a;
    }
};


class VariableDerivative : public Derivative{
private:
    int ind;

public:
    VariableDerivative(int _ind):ind(_ind){}
    
    Derivative* diffPartial(int index){
        return new ConstantDerivative(index == ind);
    }
    
    double call(VectorXd vec) const {
        return vec[ind];
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
    
    double call(VectorXd vec) const {
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
    double call(VectorXd vec) const {
        return a->call(vec) + b->call(vec);
    }
};


class DerivativeSub : public Derivative{
private:
    Derivative* a;
    Derivative* b;

public:
    DerivativeSub(Derivative* _a, Derivative* _b){a = _a, b = _b;}
    Derivative* diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) - b->call(vec);
    }
};


class DerivativeMultiply : public Derivative{
private:
    Derivative* a;
    Derivative* b;

public:
    DerivativeMultiply(Derivative* _a, Derivative* _b){a = _a, b = _b;}
    Derivative* diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) * b->call(vec);
    }
};


class DerivativeDivide : public Derivative{
private:
    Derivative* a;
    Derivative* b;

public:
    DerivativeDivide(Derivative* _a, Derivative* _b){a = _a, b = _b;}
    Derivative* diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) / b->call(vec);
    }
};


class Wrapper{
public:
    Derivative* inst;
    Wrapper():inst(nullptr){}
    Wrapper(Derivative* _inst):inst(_inst){}
    Wrapper(double x):inst(new ConstantDerivative(x)){}

    Wrapper diffPartial(int index){
        assert(inst);
        return inst->diffPartial(index);
    }

    double operator()(VectorXd vec) const {
        assert(inst);
        return inst->call(vec);
    }
};


Wrapper operator+(Wrapper a, Wrapper b){
    return new DerivativeAdd(a.inst, b.inst);
}


Wrapper operator-(Wrapper a, Wrapper b){
    return new DerivativeSub(a.inst, b.inst);
}


Wrapper operator*(Wrapper a, Wrapper b){
    return new DerivativeMultiply(a.inst, b.inst);
}


Wrapper operator/(Wrapper a, Wrapper b){
    return new DerivativeDivide(a.inst, b.inst);
}

} // Namespace Eigen
