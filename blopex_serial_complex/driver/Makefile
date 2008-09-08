include ../Makefile.inc

CPPFLAGS = -I../multivector -I../matmultivec -I../pcg_multi -I$(LOBPCG_ROOT_DIR)/include

headers = ../multivector/multi_vector.h\
../pcg_multi/pcg_multi.h ../matmultivec/matmultivec.h

objfiles = ../multivector/multi_vector.o\
../matmultivec/matmultivec.o ../pcg_multi/pcg_multi.o

lobpcglibLDFLAGS = -L$(LOBPCG_ROOT_DIR)/lib -lBLOPEX

serial_driver: serial_driver.o $(objfiles)
	$(CC) -o $@ $^ $(LAPACKLIB) $(lobpcglibLDFLAGS)
        
serial_driver.o: serial_driver.c $(headers)

clean: FORCE
	rm -f serial_driver serial_driver.o

FORCE:
