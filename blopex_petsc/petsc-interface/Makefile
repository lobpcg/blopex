include ../Makefile.inc

CFLAGS           =
CPPFLAGS         = -I$(LOBPCG_ROOT_DIR)/multivector
LOCDIR           = petsc-interface/

include ${PETSC_DIR}/bmake/common/variables

petsc-interface.o: petsc-interface.c petsc-interface.h
	$(PETSC_COMPILE_SINGLE) petsc-interface.c
        
clean: FORCE
	rm -f petsc-interface.o

FORCE:
