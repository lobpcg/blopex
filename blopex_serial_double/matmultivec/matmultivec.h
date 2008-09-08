/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ University of Colorado at Denver                                  */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#ifndef MATMULTIVECTOR_FUNCTION_PROTOTYPES
#define MATMULTIVECTOR_FUNCTION_PROTOTYPES

#include "multi_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

int serial_Lapl1DMatMultiVec( serial_Multi_Vector * x,
                            serial_Multi_Vector * y  );
                            
                            
#ifdef __cplusplus
}
#endif

#endif /*MATMULTIVECTOR_FUNCTION_PROTOTYPES */                            
