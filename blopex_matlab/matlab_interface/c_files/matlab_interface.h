/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.2) MIT / Apache-2.0 License                   */
/* @@@ Copyright 2010-2019 Team https://github.com/lobpcg/blopex       */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#ifndef MATLAB_INTERFACE_HEADER
#define MATLAB_INTERFACE_HEADER

#include "interpreter.h"
#include "fortran_matrix.h"   /* add def of komplex data type */

int MATLABSetupInterpreter( mv_InterfaceInterpreter *i );
void MatMultiVec (void * data, void * x, void * y);
BlopexInt dpotrf_redirect (char *uplo, BlopexInt *n, double *R, BlopexInt *lda, BlopexInt *info);
BlopexInt dsygv_redirect (BlopexInt *itype, char *jobz, char *uplo, BlopexInt * n,
                    double *a, BlopexInt *lda, double *b, BlopexInt *ldb,
                    double *w, double *work, BlopexInt *lwork, BlopexInt *info);
BlopexInt zpotrf_redirect (char *uplo, BlopexInt *n, komplex *R, BlopexInt *lda, BlopexInt *info);
BlopexInt zhegv_redirect (BlopexInt *itype, char *jobz, char *uplo, BlopexInt *
                    n, komplex *a, BlopexInt *lda, komplex *b, BlopexInt *ldb,
                    double *w, komplex *work, BlopexInt *lwork, double* rwork, BlopexInt *info);
#endif /* MATLAB_INTERFACE_HEADER */
