all: examples

examples: examples/many_variables.out examples/basic-calculation.out examples/interior-point-method.out examples/gauss-newton.out examples/levenberg-marquardt.out

%.out: %.cpp Derivative.cpp Derivative.h
	g++ $< -o $@ -I eigen/ -I . -std=c++11
