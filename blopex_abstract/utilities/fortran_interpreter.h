/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#ifndef LOBPCG_FORTRAN_INTERPRETER_HEADER
#define LOBPCG_FORTRAN_INTERPRETER_HEADER

/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk, Evgueni Ovtchinnikov */
/* and Don McCuan  */

typedef struct
{
  utilities_FortranMatrix* (*FortranMatrixCreate)  ( void );
  void (*FortranMatrixAllocateData) ( long, long, utilities_FortranMatrix* );
  void (*FortranMatrixWrap) ( void* , long, long, long, utilities_FortranMatrix* );
  void (*FortranMatrixDestroy) ( utilities_FortranMatrix* );
  long (*FortranMatrixGlobalHeight) ( utilities_FortranMatrix* );
  long (*FortranMatrixHeight) ( utilities_FortranMatrix* );
  long (*FortranMatrixWidth) ( utilities_FortranMatrix* );
  void* (*FortranMatrixValues) ( utilities_FortranMatrix* );
  void (*FortranMatrixClear) ( utilities_FortranMatrix* );
  void (*FortranMatrixClearL) ( utilities_FortranMatrix* );
  void (*FortranMatrixSetToIdentity) ( utilities_FortranMatrix* );
  void (*FortranMatrixTransposeSquare) ( utilities_FortranMatrix* );
  void (*FortranMatrixSymmetrize) ( utilities_FortranMatrix* );
  void (*FortranMatrixCopy) ( utilities_FortranMatrix* src, int,
                              utilities_FortranMatrix* dest);
  void (*FortranMatrixIndexCopy) ( int* index,
                                   utilities_FortranMatrix* src, int,
                                   utilities_FortranMatrix* dest);
  void (*FortranMatrixSetDiagonal) ( utilities_FortranMatrix* mtx ,
                                     utilities_FortranMatrix* d );
  void (*FortranMatrixGetDiagonal) ( utilities_FortranMatrix* mtx,
                                     utilities_FortranMatrix* d );
  void (*FortranMatrixAdd) ( double a, utilities_FortranMatrix* mtxA,
                                       utilities_FortranMatrix* mtxB,
                                       utilities_FortranMatrix* mtxC);
  void (*FortranMatrixDMultiply) ( utilities_FortranMatrix* d,
                                   utilities_FortranMatrix* mtx);
  void (*FortranMatrixMultiplyD) ( utilities_FortranMatrix* mtx,
                                   utilities_FortranMatrix* d);
  void (*FortranMatrixMultiply) ( utilities_FortranMatrix* mtxA, int tA,
                                  utilities_FortranMatrix* mtxB, int tB,
                                  utilities_FortranMatrix* mtxC );
  double (*FortranMatrixFNorm) ( utilities_FortranMatrix* mtx );
  double (*FortranMatrixAbs) ( utilities_FortranMatrix* mtx, long i, long j );

  void* (*FortranMatrixValuePtr) ( utilities_FortranMatrix* mtx, long i, long j );
  double (*FortranMatrixMaxValue) ( utilities_FortranMatrix* mtx );
  void (*FortranMatrixSelectBlock) ( utilities_FortranMatrix* mtx,
                                     long iFrom, long iTo, long jFrom, long jTo,
                                     utilities_FortranMatrix* block);
  void (*FortranMatrixUpperInv) ( utilities_FortranMatrix* u );
  int (*FortranMatrixPrint) ( utilities_FortranMatrix* mtx, char fileName[] );

} utilities_FortranInterpreter;

#endif /* LOBPCG_FORTRAN_INTERPRETER_HEADER */
