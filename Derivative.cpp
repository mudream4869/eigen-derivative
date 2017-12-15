#include "Derivative.h"

Derivative DerivativeAdd::diffPartial(int index){
    return SingleDerivative(
        [this](VectorXd vec){
            return a(vec) + b(vec);
        },
        {
            new DerivativeAdd(a.diffPartial(index), b.diffPartial(index))
        }
    );
}
