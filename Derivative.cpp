#include <functional>
#include <limits>
#include <cmath>
#include <cassert>
#include "Derivative.h"

using std::function;

namespace Eigen{

// newXXXNode implement some reduce when creating the node.

ptrDerivativeNode newDerivativeAddNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    // Very simple reduce for one of node is 0
    if(b->isConstant(0)) return a;
    if(a->isConstant(0)) return b;

    return ptrDerivativeNode(new DerivativeAddNode(a, b));
}

ptrDerivativeNode newDerivativeSubNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    // Very simple reduce for one of node is 0
    if(b->isConstant(0)) return a;

    return ptrDerivativeNode(new DerivativeSubNode(a, b));
}

ptrDerivativeNode newDerivativeMultiplyNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    // Very simple reduce for one of node is 0 or 1
    if(a->isConstant(0) or b->isConstant(0))
        return ptrDerivativeNode(new ConstantDerivativeNode(0));
    if(a->isConstant(1)) return b;
    if(b->isConstant(1)) return a;

    return ptrDerivativeNode(new DerivativeMultiplyNode(a, b));
}

ptrDerivativeNode newDerivativeDivideNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    // Very simple reduce for one of node is 0 or 1
    if(a->isConstant(0))
        return ptrDerivativeNode(new ConstantDerivativeNode(0));
    if(b->isConstant(1)) return a;

    return ptrDerivativeNode(new DerivativeDivideNode(a, b));
}

ptrDerivativeNode newDerivativePowNode(const ptrDerivativeNode& a, double p){
    return ptrDerivativeNode(new DerivativePowNode(a, p));
}

ptrDerivativeNode newDerivativeExpNode(const ptrDerivativeNode& a){
    return ptrDerivativeNode(new DerivativeExpNode(a));
}

ptrDerivativeNode newDerivativeLogNode(const ptrDerivativeNode& a){
    return ptrDerivativeNode(new DerivativeLogNode(a));
}
   

ptrDerivativeNode DerivativeNode::diffPartial(int index){
    // Save the calculated partial differential node to save time. 
    if(not dp_map.count(index))
        dp_map[index] = this->_diffPartial(index);
    return dp_map[index];
}

ptrDerivativeNode DerivativeNode::_diffPartial(int index){
    assert(0 and "DerivativeNode doesn't implement partial differential function.");
}   

double DerivativeNode::call(const VectorXd& vec) const {
    assert(0 and "DerivativeNode doesn't implement call function.");
}

void DerivativeNode::print(std::ostream& stream) const {
    assert(0 and "DerivativeNode doesn't implement print function.");
}

bool DerivativeNode::isConstant(double c) const {
    return false;
}


ConstantDerivativeNode::ConstantDerivativeNode(double _a):a(_a){}

ptrDerivativeNode ConstantDerivativeNode::_diffPartial(int index){
    return ptrDerivativeNode(new ConstantDerivativeNode(0));
}

double ConstantDerivativeNode::call(const VectorXd& vec) const {
    return a;
}

void ConstantDerivativeNode::print(std::ostream& stream) const {
    stream << a;
    return;
}

bool ConstantDerivativeNode::isConstant(double c) const {
    return std::abs(a-c) < std::numeric_limits<double>::min();
}

VariableDerivativeNode::VariableDerivativeNode(int _ind):ind(_ind){
}

ptrDerivativeNode VariableDerivativeNode::_diffPartial(int index){
    return ptrDerivativeNode(new ConstantDerivativeNode(index == ind));
}

double VariableDerivativeNode::call(const VectorXd& vec) const {
    return vec[ind];
}

void VariableDerivativeNode::print(std::ostream& stream) const {
    stream << "x[" << ind << "]";
    return;
}


LinearDerivativeNode::LinearDerivativeNode(VectorXd _v):v(_v){
}

ptrDerivativeNode LinearDerivativeNode::_diffPartial(int index){
    return ptrDerivativeNode(new ConstantDerivativeNode(v[index]));
}

double LinearDerivativeNode::call(const VectorXd& vec) const {
    return v.transpose()*vec;
}


DerivativeAddNode::DerivativeAddNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b):a(_a), b(_b){
}

double DerivativeAddNode::call(const VectorXd& vec) const {
    return a->call(vec) + b->call(vec);
}


DerivativeSubNode::DerivativeSubNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b):a(_a), b(_b){
}

double DerivativeSubNode::call(const VectorXd& vec) const {
    return a->call(vec) - b->call(vec);
}


DerivativeMultiplyNode::DerivativeMultiplyNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b):a(_a), b(_b){
}

double DerivativeMultiplyNode::call(const VectorXd& vec) const {
    return a->call(vec) * b->call(vec);
}


DerivativeDivideNode::DerivativeDivideNode(const ptrDerivativeNode& _a, const ptrDerivativeNode& _b):a(_a), b(_b){
}

double DerivativeDivideNode::call(const VectorXd& vec) const {
    return a->call(vec) / b->call(vec);
}


DerivativePowNode::DerivativePowNode(const ptrDerivativeNode& _a, double _p):a(_a), p(_p){
}

double DerivativePowNode::call(const VectorXd& vec) const {
    return std::pow(a->call(vec), p);
}


DerivativeExpNode::DerivativeExpNode(const ptrDerivativeNode& _a):a(_a){
}

double DerivativeExpNode::call(const VectorXd& vec) const {
    return std::exp(a->call(vec));
}


DerivativeLogNode::DerivativeLogNode(const ptrDerivativeNode& _a):a(_a){
}

double DerivativeLogNode::call(const VectorXd& vec) const {
    return std::log(a->call(vec));
}


Derivative::Derivative(ptrDerivativeNode _inst):inst(_inst){
}

Derivative::Derivative(double x):inst(new ConstantDerivativeNode(x)){
}

Derivative Derivative::Variable(int ind){
    return ptrDerivativeNode(new VariableDerivativeNode(ind));
}

Derivative Derivative::diffPartial(int index){
    // inst might be null
    assert(inst);
    return inst->diffPartial(index);
}

double Derivative::operator()(const VectorXd& vec) const {
    // inst might be null
    assert(inst);
    return inst->call(vec);
}


std::ostream& operator<< (std::ostream& stream, const Derivative& a){
    a.inst->print(stream);
    return stream;
}


void LinearDerivativeNode::print(std::ostream& stream) const {
    // TODO(Mudream): output more simple formula
    stream << "(";
    for(int lx = 0;lx < v.size();lx++){
        stream << v[lx] << "x[" << lx << "]";
        if(lx+1 < v.size())
            stream << " + ";
    }
    stream << ")";
    return;
}

void DerivativeAddNode::print(std::ostream& stream) const {
    stream << "("; 
    a->print(stream);
    stream << " + ";
    b->print(stream);
    stream << ")"; 
    return;
}

void DerivativeSubNode::print(std::ostream& stream) const {
    stream << "("; 
    a->print(stream);
    stream << " - ";
    b->print(stream);
    stream << ")"; 
    return;
}

void DerivativeMultiplyNode::print(std::ostream& stream) const {
    stream << "("; 
    a->print(stream);
    stream << " * ";
    b->print(stream);
    stream << ")"; 
    return;
}

void DerivativeDivideNode::print(std::ostream& stream) const {
    stream << "("; 
    a->print(stream);
    stream << " / ";
    b->print(stream);
    stream << ")"; 
    return;
}

void DerivativePowNode::print(std::ostream& stream) const {
    stream << "("; 
    a->print(stream);
    stream << "**" << p << ")";
    return;
}

void DerivativeExpNode::print(std::ostream& stream) const {
    stream << "Exp("; 
    a->print(stream);
    stream << ")"; 
    return;
}

void DerivativeLogNode::print(std::ostream& stream) const {
    stream << "Log("; 
    a->print(stream);
    stream << ")"; 
    return;
}


// Calculate the differential acording to differential rule

ptrDerivativeNode DerivativeAddNode::_diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeAddNode(ad, bd); 
}

ptrDerivativeNode DerivativeSubNode::_diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeSubNode(ad, bd);
}

ptrDerivativeNode DerivativeMultiplyNode::_diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeAddNode(
        newDerivativeMultiplyNode(ad, b),
        newDerivativeMultiplyNode(bd, a)
    );
}

ptrDerivativeNode DerivativeDivideNode::_diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeDivideNode(
        newDerivativeSubNode(
            newDerivativeMultiplyNode(ad, b),
            newDerivativeMultiplyNode(bd, a)
        ),
        newDerivativeMultiplyNode(b, b)
    );
}

ptrDerivativeNode DerivativePowNode::_diffPartial(int index){
    return newDerivativeMultiplyNode(
        ptrDerivativeNode(new ConstantDerivativeNode(p)),
        newDerivativePowNode(a, p-1)
    );
}

ptrDerivativeNode DerivativeExpNode::_diffPartial(int index){
    return newDerivativeMultiplyNode(
        a->diffPartial(index),
        newDerivativeExpNode(a)
    );
}

ptrDerivativeNode DerivativeLogNode::_diffPartial(int index){
    return newDerivativeDivideNode(
        a->diffPartial(index), a
    );
}


// Operator on Wrapper

Derivative operator+(const Derivative& a, const Derivative& b){
    return newDerivativeAddNode(a.inst, b.inst);
}

Derivative operator-(const Derivative& a, const Derivative& b){
    return newDerivativeSubNode(a.inst, b.inst);
}

Derivative operator*(const Derivative& a, const Derivative& b){
    return newDerivativeMultiplyNode(a.inst, b.inst);
}

Derivative operator/(const Derivative& a, const Derivative& b){
    return newDerivativeDivideNode(a.inst, b.inst);
}

Derivative pow(const Derivative& a, double p){
    return newDerivativePowNode(a.inst, p);
}

Derivative exp(const Derivative& a){
    return newDerivativeExpNode(a.inst);
}

Derivative log(const Derivative& a){
    return newDerivativeLogNode(a.inst);
}

} // namespace Eigen
