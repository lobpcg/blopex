#ifndef MATLAB_INTERFACE_HEADER
#define MATLAB_INTERFACE_HEADER

#include "interpreter.h"

int MATLABSetupInterpreter( mv_InterfaceInterpreter *i );
void MatMultiVec (void * data, void * x, void * y);

#endif /* MATLAB_INTERFACE_HEADER */
