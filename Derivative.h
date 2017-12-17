#include <Eigen/Dense>
#include <functional>
#include <vector>

using Eigen::VectorXd;
using std::function;
using std::vector;

typedef function<double(VectorXd)> MultiVar;

namespace Eigen{

class DerivativeNode{
public:
    DerivativeNode(){}
    virtual DerivativeNode* diffPartial(int index){}
    virtual double call(VectorXd vec) const {}
};


class ConstantDerivativeNode : public DerivativeNode{
private:
    double a;

public:
    ConstantDerivativeNode(double _a):a(_a){}
    
    DerivativeNode* diffPartial(int index){
        return new ConstantDerivativeNode(0);
    }
    
    double call(VectorXd vec) const {
        return a;
    }
};


class VariableDerivativeNode : public DerivativeNode{
private:
    int ind;

public:
    VariableDerivativeNode(int _ind):ind(_ind){}
    
    DerivativeNode* diffPartial(int index){
        return new ConstantDerivativeNode(index == ind);
    }
    
    double call(VectorXd vec) const {
        return vec[ind];
    }
};



class LinearDerivativeNode : public DerivativeNode{
private:
    VectorXd v;

public:
    LinearDerivativeNode(VectorXd _v){
        v = _v;
    }
 
    DerivativeNode* diffPartial(int index){
        return new ConstantDerivativeNode(v[index]);
    }
    
    double call(VectorXd vec) const {
        return v.transpose()*vec;
    }
};


class DerivativeAddNode : public DerivativeNode{
private:
    DerivativeNode* a;
    DerivativeNode* b;

public:
    DerivativeAddNode(DerivativeNode* _a, DerivativeNode* _b){a = _a, b = _b;}
    DerivativeNode* diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) + b->call(vec);
    }
};


class DerivativeSubNode : public DerivativeNode{
private:
    DerivativeNode* a;
    DerivativeNode* b;

public:
    DerivativeSubNode(DerivativeNode* _a, DerivativeNode* _b){a = _a, b = _b;}
    DerivativeNode* diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) - b->call(vec);
    }
};


class DerivativeMultiplyNode : public DerivativeNode{
private:
    DerivativeNode* a;
    DerivativeNode* b;

public:
    DerivativeMultiplyNode(DerivativeNode* _a, DerivativeNode* _b){a = _a, b = _b;}
    DerivativeNode* diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) * b->call(vec);
    }
};


class DerivativeDivideNode : public DerivativeNode{
private:
    DerivativeNode* a;
    DerivativeNode* b;

public:
    DerivativeDivideNode(DerivativeNode* _a, DerivativeNode* _b){a = _a, b = _b;}
    DerivativeNode* diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) / b->call(vec);
    }
};


class Wrapper{
public:
    DerivativeNode* inst;
    Wrapper(DerivativeNode* _inst = nullptr):inst(_inst){}
    Wrapper(double x):inst(new ConstantDerivativeNode(x)){}

    static Wrapper Variable(int ind){
        return new VariableDerivativeNode(ind);
    }

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
    return new DerivativeAddNode(a.inst, b.inst);
}


Wrapper operator-(Wrapper a, Wrapper b){
    return new DerivativeSubNode(a.inst, b.inst);
}


Wrapper operator*(Wrapper a, Wrapper b){
    return new DerivativeMultiplyNode(a.inst, b.inst);
}


Wrapper operator/(Wrapper a, Wrapper b){
    return new DerivativeDivideNode(a.inst, b.inst);
}

} // Namespace Eigen
