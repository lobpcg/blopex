/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ University of Colorado at Denver                                  */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#ifndef HYPRE_PARCSR_INTERFACE_INTERPRETER
#define HYPRE_PARCSR_INTERFACE_INTERPRETER

#include "interpreter.h"
#include "HYPRE_MatvecFunctions.h"

#ifdef __cplusplus
extern "C" {
#endif

int
HYPRE_SStructSetupInterpreter( mv_InterfaceInterpreter *i );

int
HYPRE_SStructSetupMatvec(HYPRE_MatvecFunctions * mv);

#ifdef __cplusplus
}
#endif

#endif
