#include "Derivative.h"

namespace Eigen{

DerivativeNode* newDerivativeAddNode(DerivativeNode* a, DerivativeNode* b){
    if(b->isConstant(0)) return a;
    if(a->isConstant(0)) return b;

    return new DerivativeAddNode(a, b);
}

DerivativeNode* newDerivativeSubNode(DerivativeNode* a, DerivativeNode* b){
    if(b->isConstant(0)) return a;

    return new DerivativeSubNode(a, b);
}

DerivativeNode* newDerivativeMultiplyNode(DerivativeNode* a, DerivativeNode* b){
    if(a->isConstant(0) or b->isConstant(0))
        return new ConstantDerivativeNode(0);
    if(a->isConstant(1)) return b;
    if(b->isConstant(1)) return a;

    return new DerivativeMultiplyNode(a, b);
}

DerivativeNode* newDerivativeDivideNode(DerivativeNode* a, DerivativeNode* b){
    if(a->isConstant(0))
        return new ConstantDerivativeNode(0);
    if(b->isConstant(1)) return a;

    return new DerivativeDivideNode(a, b);
}


DerivativeNode* DerivativeAddNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeAddNode(ad, bd); 
}

DerivativeNode* DerivativeSubNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeSubNode(ad, bd);
}

DerivativeNode* DerivativeMultiplyNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeAddNode(
        newDerivativeMultiplyNode(ad, b),
        newDerivativeMultiplyNode(bd, a)
    );
}

DerivativeNode* DerivativeDivideNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return newDerivativeDivideNode(
        newDerivativeSubNode(
            newDerivativeMultiplyNode(ad, b),
            newDerivativeMultiplyNode(bd, a)
        ),
        newDerivativeMultiplyNode(b, b)
    );
}

Derivative operator+(Derivative a, Derivative b){
    return newDerivativeAddNode(a.inst, b.inst);
}

Derivative operator-(Derivative a, Derivative b){
    return newDerivativeSubNode(a.inst, b.inst);
}

Derivative operator*(Derivative a, Derivative b){
    return newDerivativeMultiplyNode(a.inst, b.inst);
}

Derivative operator/(Derivative a, Derivative b){
    return newDerivativeDivideNode(a.inst, b.inst);
}

} // Namespace Eigen
