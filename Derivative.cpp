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


ptrDerivativeNode DerivativeAddNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeAddNode(ad, bd); 
}

ptrDerivativeNode DerivativeSubNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeSubNode(ad, bd);
}

ptrDerivativeNode DerivativeMultiplyNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeAddNode(
        newDerivativeMultiplyNode(ad, b),
        newDerivativeMultiplyNode(bd, a)
    );
}

ptrDerivativeNode DerivativeDivideNode::diffPartial(int index){
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
