#include "Derivative.h"

namespace Eigen{

ptrDerivativeNode newDerivativeAddNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    if(b->isConstant(0)) return a;
    if(a->isConstant(0)) return b;

    return ptrDerivativeNode(new DerivativeAddNode(a, b));
}

ptrDerivativeNode newDerivativeSubNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    if(b->isConstant(0)) return a;

    return ptrDerivativeNode(new DerivativeSubNode(a, b));
}

ptrDerivativeNode newDerivativeMultiplyNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    if(a->isConstant(0) or b->isConstant(0))
        return ptrDerivativeNode(new ConstantDerivativeNode(0));
    if(a->isConstant(1)) return b;
    if(b->isConstant(1)) return a;

    return ptrDerivativeNode(new DerivativeMultiplyNode(a, b));
}

ptrDerivativeNode newDerivativeDivideNode(const ptrDerivativeNode& a, const ptrDerivativeNode& b){
    if(a->isConstant(0))
        return ptrDerivativeNode(new ConstantDerivativeNode(0));
    if(b->isConstant(1)) return a;

    return ptrDerivativeNode(new DerivativeDivideNode(a, b));
}

void LinearDerivativeNode::print(std::ostream& stream) const {
    stream << "(";
    for(int lx = 0;lx < v.size();lx++){
        stream << v[lx] << "x[" << lx << "]";
        if(lx+1 < v.size())
            stream << " + ";
    }
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

} // Namespace Eigen
