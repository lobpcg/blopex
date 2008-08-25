VPATH = driver:matlab_interface/c_files
CC = mex
CFLAGS = -g
CPPFLAGS = -Imatlab_interface/c_files -I../blopex_abstract/include
LIB_DIRS = -L../blopex_abstract/lib
LIBS =  -lBLOPEX

driver/blopex_matlab_gateway.mexglx: blopex_matlab_gateway.c multivector_for_matlab.h matlab_interface.h \
                                     multivector_for_matlab.o matlab_interface.o
	$(CC) $(CFLAGS) $(CPPFLAGS) $< matlab_interface/c_files/multivector_for_matlab.o \
         matlab_interface/c_files/matlab_interface.o $(LIB_DIRS) $(LIBS) -outdir $(@D)

matlab_interface/c_files/multivector_for_matlab.o: multivector_for_matlab.c multivector_for_matlab.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -outdir $(@D)

matlab_interface/c_files/matlab_interface.o: matlab_interface.c matlab_interface.h multivector_for_matlab.h
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -outdir $(@D)

clean: FORCE
	rm -f matlab_interface/c_files/multivector_for_matlab.o matlab_interface/c_files/matlab_interface.o \
              driver/blopex_matlab_gateway.mexglx

FORCE:
