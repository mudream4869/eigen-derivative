#include <Eigen/Dense>
#include <functional>
#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> 

using Eigen::VectorXd;
using std::function;
using std::vector;

typedef function<double(VectorXd)> MultiVar;

namespace Eigen{

class DerivativeNode;
typedef std::shared_ptr<DerivativeNode> ptrDerivativeNode;

class DerivativeNode{
public:
    DerivativeNode(){}
    virtual ptrDerivativeNode diffPartial(int index){
        assert(0 and "DerivativeNode doesn't implement.");
    }

    virtual double call(VectorXd vec) const {
        assert(0 and "DerivativeNode doesn't implement.");
    }

    virtual void print(std::ostream& stream) const {
        assert(0 and "DerivativeNode doesn't implement.");
    }

    virtual bool isConstant(double c) const { return false; }
};


class ConstantDerivativeNode : public DerivativeNode{
private:
    double a;

public:
    ConstantDerivativeNode(double _a):a(_a){}
    
    ptrDerivativeNode diffPartial(int index){
        return ptrDerivativeNode(new ConstantDerivativeNode(0));
    }
    
    double call(VectorXd vec) const {
        return a;
    }

    void print(std::ostream& stream) const {
        stream << a;
        return;
    }

    bool isConstant(double c) const {
        return std::abs(a-c) < std::numeric_limits<double>::min();
    }
};


class VariableDerivativeNode : public DerivativeNode{
private:
    int ind;

public:
    VariableDerivativeNode(int _ind):ind(_ind){}
    
    ptrDerivativeNode diffPartial(int index){
        return ptrDerivativeNode(new ConstantDerivativeNode(index == ind));
    }
    
    double call(VectorXd vec) const {
        return vec[ind];
    }

    void print(std::ostream& stream) const {
        stream << "x[" << ind << "]";
        return;
    }
};



class LinearDerivativeNode : public DerivativeNode{
private:
    VectorXd v;

public:
    LinearDerivativeNode(VectorXd _v){
        v = _v;
    }
 
    ptrDerivativeNode diffPartial(int index){
        return ptrDerivativeNode(new ConstantDerivativeNode(v[index]));
    }
    
    double call(VectorXd vec) const {
        return v.transpose()*vec;
    }

    void print(std::ostream& stream) const;
};


class DerivativeAddNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeAddNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b){a = _a, b = _b;}
    ptrDerivativeNode diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) + b->call(vec);
    }

    void print(std::ostream& stream) const;
};


class DerivativeSubNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeSubNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b){a = _a, b = _b;}
    ptrDerivativeNode diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) - b->call(vec);
    }

    void print(std::ostream& stream) const;
};


class DerivativeMultiplyNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeMultiplyNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b){a = _a, b = _b;}
    ptrDerivativeNode diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) * b->call(vec);
    }

    void print(std::ostream& stream) const;
};


class DerivativeDivideNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeDivideNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b){a = _a, b = _b;}
    ptrDerivativeNode diffPartial(int index);
    double call(VectorXd vec) const {
        return a->call(vec) / b->call(vec);
    }

    void print(std::ostream& stream) const;
};


class Derivative{
public:
    ptrDerivativeNode inst;
    Derivative(ptrDerivativeNode _inst = nullptr):inst(_inst){}
    Derivative(double x):inst(new ConstantDerivativeNode(x)){}

    static Derivative Variable(int ind){
        return ptrDerivativeNode(new VariableDerivativeNode(ind));
    }

    Derivative diffPartial(int index){
        assert(inst);
        return inst->diffPartial(index);
    }

    double operator()(VectorXd vec) const {
        assert(inst);
        return inst->call(vec);
    }
};

std::ostream& operator<< (std::ostream& stream, const Derivative& a){
    a.inst->print(stream);
    return stream;
}

Derivative operator+(const Derivative& a, const Derivative& b);
Derivative operator-(const Derivative& a, const Derivative& b);
Derivative operator*(const Derivative& a, const Derivative& b);
Derivative operator/(const Derivative& a, const Derivative& b);

} // Namespace Eigen
