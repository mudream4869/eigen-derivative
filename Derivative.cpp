#include "Derivative.h"

Derivative* DerivativeAdd::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new SingleDerivative(
        [ad, bd](VectorXd vec){
            return ad->call(vec) + bd->call(vec);
        },
        [ad, bd](int i){
            return new DerivativeAdd(ad->diffPartial(i), bd->diffPartial(i));
        }
    );
}

Derivative* DerivativeSub::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new SingleDerivative(
        [ad, bd](VectorXd vec){
            return ad->call(vec) - bd->call(vec);
        },
        [ad, bd](int i){
            return new DerivativeSub(ad->diffPartial(i), bd->diffPartial(i));
        }
    );
}

Derivative* DerivativeMultiply::diffPartial(int index){
    auto ad = a->diffPartial(index), bd = b->diffPartial(index);
    return new SingleDerivative(
        [ad, bd](VectorXd vec){
            return ad->call(vec) * bd->call(vec);
        },
        [ad, bd](int i){
            return new DerivativeAdd(
                new DerivativeMultiply(ad->diffPartial(i), bd),
                new DerivativeMultiply(bd->diffPartial(i), ad)
            );
        }
    );
}
