all: tests examples

tests: tests/many_variables.out tests/local_test.out tests/eigen.out

examples: examples/basic-calculation.out examples/interior-point-method.out examples/gauss-newton.out examples/levenberg-marquardt.out

%.out: %.cpp Derivative.o Derivative.h
	g++ Derivative.o $< -o $@ -I eigen/ -I . -std=c++11

Derivative.o: Derivative.cpp
	g++ Derivative.cpp -I eigen/ -I . -std=c++11 -c

clear:
	rm tests/*.out examples/*.out
