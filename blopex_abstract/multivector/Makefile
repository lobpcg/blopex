include ../Makefile.inc

all: multivector.o temp_multivector.o

multivector.o: multivector.c multivector.h
temp_multivector.o: temp_multivector.c temp_multivector.h interpreter.h


clean: FORCE
	rm -f multivector.o temp_multivector.o

FORCE:
