/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#ifndef PETSC_INTERFACE_HEADER
#define PETSC_INTERFACE_HEADER

#include "interpreter.h"

int PETSC_dpotrf_interface (char *uplo, int *n, double *a, int * lda, int *info);

int PETSC_dsygv_interface (int *itype, char *jobz, char *uplo, int *
                    n, double *a, int *lda, double *b, int *ldb,
                    double *w, double *work, int *lwork, int *info);

int PETSC_zpotrf_interface (char *uplo, int *n, komplex *a, int * lda, int *info);

int PETSC_zsygv_interface (int *itype, char *jobz, char *uplo, int *
                    n, komplex *a, int *lda, komplex *b, int *ldb,
                    double *w, komplex *work, int *lwork, double *rwork, int *info);

void *
PETSC_MimicVector( void *vvector );

int
PETSC_DestroyVector( void *vvector );

int
PETSC_InnerProd( void *x, void *y, void *result );

int
PETSC_CopyVector( void *x, void *y );

int
PETSC_ClearVector( void *x );

int
PETSC_SetRandomValues( void* v, int seed );

int
PETSC_ScaleVector( void *alpha, void   *x);

int
PETSC_Axpy( void *alpha,
                void   *x,
                void   *y );

int
LOBPCG_InitRandomContext(void);

int
LOBPCG_DestroyRandomContext(void);

int
PETSCSetupInterpreter( mv_InterfaceInterpreter *ii );

#endif /* PETSC_INTERFACE_HEADER */
