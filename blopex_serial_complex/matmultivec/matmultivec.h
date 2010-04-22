/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1) LGPL Version 2.1 or above.See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team http://code.google.com/p/blopex/     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#ifndef BlopexInt
#define BlopexInt int
#endif

#ifndef MATMULTIVECTOR_FUNCTION_PROTOTYPES
#define MATMULTIVECTOR_FUNCTION_PROTOTYPES

#include "multi_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

BlopexInt serial_Lapl1DMatMultiVec( serial_Multi_Vector * x,
                            serial_Multi_Vector * y  );
                            
                            
#ifdef __cplusplus
}
#endif

#endif /*MATMULTIVECTOR_FUNCTION_PROTOTYPES */                            
