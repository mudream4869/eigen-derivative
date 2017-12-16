#include "Derivative.h"

Eigen::Derivative* Eigen::DerivativeAdd::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAdd(ad, bd); 
}

Eigen::Derivative* Eigen::DerivativeSub::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeSub(ad, bd);
}

Eigen::Derivative* Eigen::DerivativeMultiply::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAdd(
        new DerivativeMultiply(ad, b),
        new DerivativeMultiply(bd, a)
    );
}

Eigen::Derivative* Eigen::DerivativeDivide::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeDivide(
        new DerivativeSub(
            new DerivativeMultiply(ad, b),
            new DerivativeMultiply(bd, a)
        ),
        new DerivativeMultiply(b, b)
    );
}
