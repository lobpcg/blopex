/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */
/* and Don McCuan  */

#ifndef LOBPCG_INTERFACE_INTERPRETER
#define LOBPCG_INTERFACE_INTERPRETER

typedef struct
{
  /* vector operations */
  void*  (*CreateVector)  ( void *vector );
  int    (*DestroyVector) ( void *vector );

  int    (*InnerProd)     ( void *x, void *y, void *result );
  int    (*CopyVector)    ( void *x, void *y );
  int    (*ClearVector)   ( void *x );
  int    (*SetRandomValues)   ( void *x, int seed );
  int    (*ScaleVector)   ( double alpha, void *x );
  int    (*Axpy)          ( void * alpha, void *x, void *y );
  int    (*VectorSize)    (void * vector);

  /* multivector operations */
  /* do we need the following entry? */
  void*  (*CreateMultiVector)  ( void*, int n, void *vector );
  void*  (*CopyCreateMultiVector)  ( void *x, int );
  void    (*DestroyMultiVector) ( void *x );

  int    (*Width)  ( void *x );
  int    (*Height) ( void *x );

  void   (*SetMask) ( void *x, int *mask );

  void   (*CopyMultiVector)    ( void *x, void *y );
  void   (*ClearMultiVector)   ( void *x );
  void   (*SetRandomVectors)   ( void *x, int seed );
  void   (*MultiInnerProd)     ( void *x, void *y, int, int, int, void* );
  void   (*MultiInnerProdDiag) ( void *x, void *y, int*, int, void* );
  void   (*MultiVecMat)        ( void *x, int, int, int, void*, void *y );
  void   (*MultiVecMatDiag)    ( void *x, int*, int, void*, void *y );
  void   (*MultiAxpy)          ( double alpha, void *x, void *y );
  void   (*MultiPrint)         ( void *x, char * tag,int limit );

  /* do we need the following 2 entries? */
  void   (*MultiXapy)          ( void *x, int, int, int, void*, void *y );
  void   (*Eval)               ( void (*f)( void*, void*, void* ), void*, void *x, void *y );

} mv_InterfaceInterpreter;

#endif
