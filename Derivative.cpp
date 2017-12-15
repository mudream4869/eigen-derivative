#include "Derivative.h"

Derivative* DerivativeAdd::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAdd(ad, bd); 
}

Derivative* DerivativeSub::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeSub(ad, bd);
}

Derivative* DerivativeMultiply::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeAdd(
        new DerivativeMultiply(ad, b),
        new DerivativeMultiply(bd, a)
    );
}

Derivative* DerivativeDivide::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new DerivativeDivide(
        new DerivativeSub(
            new DerivativeMultiply(ad, b),
            new DerivativeMultiply(bd, a)
        ),
        new DerivativeMultiply(b, b)
    );
}
