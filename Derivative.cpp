#include "Derivative.h"

Derivative* DerivativeAdd::diffPartial(int index){
    return new SingleDerivative(
        [this](VectorXd vec){
            return a->call(vec) + b->call(vec);
        },
        {
            new DerivativeAdd(a->diffPartial(index), b->diffPartial(index))
        }
    );
}

Derivative* DerivativeMultiply::diffPartial(int index){
    return new SingleDerivative(
        [this](VectorXd vec){
            return a->call(vec) + b->call(vec);
        },
        {
            new DerivativeAdd(
                new DerivativeMultiply(a->diffPartial(index), b),
                new DerivativeMultiply(b->diffPartial(index), a)
            )
        }
    );
}
