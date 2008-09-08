include ../Makefile.inc

CPPFLAGS = -I../multivector -I../matmultivec -I$(LOBPCG_ROOT_DIR)/include

pcg_multi.o: pcg_multi.c pcg_multi.h ../matmultivec/matmultivec.h ../multivector/multi_vector.h

clean: FORCE
	rm -f pcg_multi.o

FORCE:
