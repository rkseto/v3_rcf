SHELL := /bin/bash # Use bash syntax

default: libcjam.a

libcjam.a: cjam.f cjam0.f
	gfortran -m32 -O3 -c cjam.f -o cjam.o -I../jam/src
	gfortran -m32 -O3 -c cjam0.f -o cjam0.o -I../jam/src
	ar r libcjam.a cjam.o cjam0.o

clean:
	rm -f libcjam.a cjam.o cjam0.o
