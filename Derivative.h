#ifndef DERIVATIVE_H_
#define DERIVATIVE_H_

#include <Eigen/Dense>
#include <iostream>
#include <memory> 
#include <map>

//using Eigen::VectorXd;

namespace Eigen{

class DerivativeNode;

// Use std's shared pointer
typedef std::shared_ptr<DerivativeNode> ptrDerivativeNode;


// DerivativeNode is the base class of all the class that can do partial
// differential. It would not be used directively.
//
// The class inherit the DerivativeNode should implement three function:
// 1. _diffPartial: Do partial differential.
// 2. call: As a scalar function, calculate the value and return.
// 3. print: Use ostream to output.
class DerivativeNode{
private:
    // Save the calculated partial differential node to save time. 
    std::map<int, ptrDerivativeNode> dp_map;

public:
    DerivativeNode(){}
    ptrDerivativeNode diffPartial(int index);

    virtual ptrDerivativeNode _diffPartial(int index);
    virtual double call(const VectorXd& vec) const;
    virtual void print(std::ostream& stream) const;
    virtual bool isConstant(double c) const; 
};


// Constant Node. Sample usage:
//   ptrDerivativeNode constant_node(ConstantDerivativeNode(1.2))
class ConstantDerivativeNode : public DerivativeNode{
private:
    double a;

public:
    ConstantDerivativeNode(double _a);
    
    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const; 
    void print(std::ostream& stream) const; 
    bool isConstant(double c) const;
};


// Single Variable. Sample usgae:
//   ptrDerivateNode var(VariableDerivativeNode(2));
class VariableDerivativeNode : public DerivativeNode{
private:
    // The index of variable
    int ind;

public:
    VariableDerivativeNode(int _ind);
    
    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const; 
    void print(std::ostream& stream) const;
};


// F(x) = [v1 v2 v3]*x. Sample usage:
//   Eigen::VectorXd vec(3);
//   vec << 1, 3, 4;
//   ptrDerivateNode linear_function(LinearDerivativeNode(vec));
class LinearDerivativeNode : public DerivativeNode{
private:
    VectorXd v;

public:
    LinearDerivativeNode(VectorXd _v);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


class DerivativeAddNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeAddNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


class DerivativeSubNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeSubNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


class DerivativeMultiplyNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeMultiplyNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


class DerivativeDivideNode : public DerivativeNode{
private:
    ptrDerivativeNode a, b;

public:
    DerivativeDivideNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


class DerivativePowNode : public DerivativeNode{
private:
    ptrDerivativeNode a;
    double p;

public:
    DerivativePowNode(const ptrDerivativeNode& _a, double _p);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


// Exp(f(x))
class DerivativeExpNode : public DerivativeNode{
private:
    ptrDerivativeNode a;

public:
    DerivativeExpNode(const ptrDerivativeNode& _a);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


// Take nature log on a DerivativeNode
class DerivativeLogNode : public DerivativeNode{
private:
    ptrDerivativeNode a;

public:
    DerivativeLogNode(const ptrDerivativeNode& _a);

    ptrDerivativeNode _diffPartial(int index);
    double call(const VectorXd& vec) const;
    void print(std::ostream& stream) const;
};


// DerivativeNode pointer's wrapper
class Derivative{
public:
    ptrDerivativeNode inst;
    Derivative(ptrDerivativeNode _inst = nullptr);
    Derivative(double x);
    
    // Useful in demo ... 
    static Derivative Variable(int ind);

    Derivative diffPartial(int index);
    double operator()(const VectorXd& vec) const;
};


std::ostream& operator<< (std::ostream& stream, const Derivative& a);


// Operator on Wrapper
Derivative operator+(const Derivative& a, const Derivative& b);
Derivative operator-(const Derivative& a, const Derivative& b);
Derivative operator*(const Derivative& a, const Derivative& b);
Derivative operator/(const Derivative& a, const Derivative& b);
Derivative exp(const Derivative& a);
Derivative log(const Derivative& a);
Derivative pow(const Derivative& a, double p);

} // namespace Eigen

#endif // DERIVATIVE_H_
