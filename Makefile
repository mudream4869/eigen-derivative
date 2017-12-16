all: examples

examples: basic-calculation.out

basic-calculation.out: examples/basic-calculation.cpp Derivative.cpp Derivative.h
	g++ examples/basic-calculation.cpp -o basic-calculation.out -I eigen/ -I . -std=c++11
