#include "Derivative.h"

Eigen::DerivativeNode* Eigen::DerivativeAddNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAddNode(ad, bd); 
}

Eigen::DerivativeNode* Eigen::DerivativeSubNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeSubNode(ad, bd);
}

Eigen::DerivativeNode* Eigen::DerivativeMultiplyNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAddNode(
        new DerivativeMultiplyNode(ad, b),
        new DerivativeMultiplyNode(bd, a)
    );
}

Eigen::DerivativeNode* Eigen::DerivativeDivideNode::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeDivideNode(
        new DerivativeSubNode(
            new DerivativeMultiplyNode(ad, b),
            new DerivativeMultiplyNode(bd, a)
        ),
        new DerivativeMultiplyNode(b, b)
    );
}
