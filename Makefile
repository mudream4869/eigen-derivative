all: examples

examples: basic-calculation.out interior-point-method.out

basic-calculation.out: examples/basic-calculation.cpp Derivative.cpp Derivative.h
	g++ examples/basic-calculation.cpp -o basic-calculation.out -I eigen/ -I . -std=c++11

interior-point-method.out: examples/interior-point-method.cpp Derivative.cpp Derivative.h
	g++ examples/interior-point-method.cpp -o interior-point-method.out -I eigen/ -I . -std=c++11 
