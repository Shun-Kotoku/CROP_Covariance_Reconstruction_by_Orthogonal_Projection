test.out: MC.cpp Data.h Initialization.h Memory.h
	g++ MC.cpp -O3 -o test.out -larmadillo