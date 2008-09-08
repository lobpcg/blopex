include ../Makefile.inc

#CPPFLAGS = -I../petsc-interface -I$(LOBPCG_ROOT_DIR)/include

CFLAGS = -I../petsc-interface -I$(LOBPCG_ROOT_DIR)/include

#LOCDIR = driver/

# include ${PETSC_DIR}/bmake/common/variables
include ${PETSC_DIR}/conf/variables

driver_fiedler: driver_fiedler.o
	${CLINKER} -o driver_fiedler driver_fiedler.o ../petsc-interface/petsc-interface.o -L$(LOBPCG_ROOT_DIR)/lib -lBLOPEX ${PETSC_KSP_LIB}

driver_fiedler.o: driver_fiedler.c ../petsc-interface/petsc-interface.h
	$(PETSC_COMPILE_SINGLE) driver_fiedler.c

clean: FORCE
	rm -f driver_fiedler.o driver_fiedler

FORCE:
