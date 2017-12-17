#include "Derivative.h"

namespace Eigen{

DerivativeNode* DerivativeAddNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAddNode(ad, bd); 
}

DerivativeNode* DerivativeSubNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeSubNode(ad, bd);
}

DerivativeNode* DerivativeMultiplyNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAddNode(
        new DerivativeMultiplyNode(ad, b),
        new DerivativeMultiplyNode(bd, a)
    );
}

DerivativeNode* DerivativeDivideNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeDivideNode(
        new DerivativeSubNode(
            new DerivativeMultiplyNode(ad, b),
            new DerivativeMultiplyNode(bd, a)
        ),
        new DerivativeMultiplyNode(b, b)
    );
}

Derivative operator+(Derivative a, Derivative b){
    return new DerivativeAddNode(a.inst, b.inst);
}

Derivative operator-(Derivative a, Derivative b){
    return new DerivativeSubNode(a.inst, b.inst);
}

Derivative operator*(Derivative a, Derivative b){
    return new DerivativeMultiplyNode(a.inst, b.inst);
}

Derivative operator/(Derivative a, Derivative b){
    return new DerivativeDivideNode(a.inst, b.inst);
}

} // Namespace Eigen
