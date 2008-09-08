include ../Makefile.inc

CFLAGS   = -I../petsc-interface -I$(LOBPCG_ROOT_DIR)/include 
FFLAGS   =
SOURCEC  = petsc-interface.c
SOURCEF  =
SOURCEH  = petsc-interface.h
OBJSC    = petsc-interface.o
OBJSF    =
LIBBASE  = libpetscksp
DIRS     =
MANSEC   = 
LOCDIR   = petsc-interface/

include ${PETSC_DIR}/conf/variables
# include ${PETSC_DIR}/bmake/common/variables

petsc-interface.o: petsc-interface.c petsc-interface.h
	$(PETSC_COMPILE_SINGLE) petsc-interface.c

clean: FORCE
	rm -f petsc-interface.o

FORCE:
	
