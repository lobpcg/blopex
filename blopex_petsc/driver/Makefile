include ../Makefile.inc

CPPFLAGS         = -I../petsc-interface -I$(LOBPCG_ROOT_DIR)/include
LOCDIR           = driver/

include ${PETSC_DIR}/bmake/common/variables

driver: driver.o
	${CLINKER} -o driver driver.o ../petsc-interface/petsc-interface.o \
-L$(LOBPCG_ROOT_DIR)/lib -lBLOPEX ${PETSC_KSP_LIB}

driver.o: driver.c ../petsc-interface/petsc-interface.h
	$(PETSC_COMPILE_SINGLE) driver.c

clean: FORCE
	rm -f driver.o driver

FORCE:
