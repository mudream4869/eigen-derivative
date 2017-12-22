#include "Derivative.h"

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

} // Namespace Eigen
