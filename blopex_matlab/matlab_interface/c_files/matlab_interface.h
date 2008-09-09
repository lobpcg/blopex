/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#ifndef MATLAB_INTERFACE_HEADER
#define MATLAB_INTERFACE_HEADER

#include "interpreter.h"

int MATLABSetupInterpreter( mv_InterfaceInterpreter *i );
void MatMultiVec (void * data, void * x, void * y);

#endif /* MATLAB_INTERFACE_HEADER */
