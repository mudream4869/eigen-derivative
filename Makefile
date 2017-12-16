all: examples

examples: examples/basic-calculation.out examples/interior-point-method.out examples/gauss-newton.out

%.out: %.cpp Derivative.cpp Derivative.h
	g++ $< -o $@ -I eigen/ -I . -std=c++11
