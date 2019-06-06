/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.2) MIT / Apache-2.0 License                   */
/* @@@ Copyright 2010-2019 Team https://github.com/lobpcg/blopex       */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

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
