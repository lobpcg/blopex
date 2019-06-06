/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.2) MIT / Apache-2.0 License                   */
/* @@@ Copyright 2010-2019 Team https://github.com/lobpcg/blopex       */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/* this file contains interface routines to matlab for blopex */
#include "matlab_interface.h"

#include "interpreter.h"

#include "multivector_for_matlab.h"

#ifdef BLOPEX_MATLAB_32
 int zpotrf_ (char *uplo, int *n,
                    komplex *a, int *lda, int *info);
 int zhegv_  (int *itype, char *jobz, char *uplo, int * n,
                    komplex *a, int *lda, komplex *b, int *ldb,
                    double *w, komplex *work, int *lwork,
                    double * rwork,int *info);
 int dsygv_ (int *itype, char *jobz, char *uplo, int *
                     n, double *a, int *lda, double *b, int *ldb,
                     double *w, double *work, int *lwork, int *info);

 int dpotrf_ (char *uplo, int *n, double *a, int *
                     lda, int *info);
#else
#include "lapack.h"
#endif

void*
CreateCopyMultiVector( void* src_, BlopexInt copyValues )
{
  return MatlabMultiVectorCopyCreate((matlabMultiVectorPtr) src_,
                                     copyValues);
}

void
DestroyMultiVector( void *vvector )
{
  MatlabMultiVectorDestroy((matlabMultiVectorPtr) vvector);
}


BlopexInt
MultiVectorWidth( void* v )
{
  return MatlabMultiVectorWidth((matlabMultiVectorPtr) v);
}

void
MultiSetMask( void *vector, BlopexInt *mask )
{
  MatlabMultiVectorSetMask((matlabMultiVectorPtr) vector, mask);
}


void
CopyMultiVector( void *x, void *y)
{
  MatlabMultiVectorCopy((matlabMultiVectorPtr)x, (matlabMultiVectorPtr)y);
}


void
MultiInnerProd(void * x_, void * y_,
				    BlopexInt gh, BlopexInt h, BlopexInt w, void* v )
{
  MatlabMultiVectorInnerProd((matlabMultiVectorPtr) x_, 
                             (matlabMultiVectorPtr) y_, gh, h, w, v);
}

void
MultiInnerProdDiag( void* x_, void* y_,
					BlopexInt* mask, BlopexInt n, void* diag )
{
  MatlabMultiVectorInnerProdDiag((matlabMultiVectorPtr) x_,
                                 (matlabMultiVectorPtr) y_, mask, n, diag);
}


void
MultiVectorByDiagonal( void* x, 
				      BlopexInt* mask, BlopexInt n, void* diag,
				      void* y )
{
  MatlabMultiVectorByDiag((matlabMultiVectorPtr) x, mask, n, diag, 
                          (matlabMultiVectorPtr) y);
}


void 
MultiVectorByMatrix( void* x, 
			       BlopexInt gh, BlopexInt h, BlopexInt w, void* v,
			       void* y )
{
  MatlabMultiVectorByMatrix((matlabMultiVectorPtr) x, gh, h, w, v, 
                            (matlabMultiVectorPtr) y);
}

void
MultiVectorAxpy( double alpha, void   *x, void   *y)
{
  MatlabMultiVectorAxpy(alpha, (matlabMultiVectorPtr) x, 
                        (matlabMultiVectorPtr) y);
}

void MatMultiVec (void * data, void * x, void * y)
{
   MatlabMultiVectorMatMultiVec( (mxArray*) data, (matlabMultiVectorPtr) x,
                                 (matlabMultiVectorPtr) y  );
}

/* redirect dpotrf to matlab lapack version where BlopexInt variables are
 * defined as ptrdiff_t (same as mwSignedIndex) */

BlopexInt dpotrf_redirect (char *uplo, BlopexInt *n, double *R, BlopexInt *lda, BlopexInt *info)
{
#ifdef BLOPEX_MATLAB_32
  int xn, xlda, xinfo;
#else
  mwSignedIndex xn, xlda, xinfo;
#endif
  xn = *n;
  xlda = *lda;
  dpotrf_(uplo, &xn, R, &xlda, &xinfo);
  *info = (BlopexInt)xinfo;
}

BlopexInt dsygv_redirect (BlopexInt *itype, char *jobz, char *uplo, BlopexInt *
                    n, double *a, BlopexInt *lda, double *b, BlopexInt *ldb,
                    double *w, double *work, BlopexInt *lwork, BlopexInt *info)
{
#ifdef BLOPEX_MATLAB_32
  int xitype, xn, xlda, xldb, xlwork, xinfo;
#else
  mwSignedIndex xitype, xn, xlda, xldb, xlwork, xinfo;
#endif
  xitype = *itype;
  xn = *n;
  xlda = *lda;
  xldb = *ldb;
  xlwork = *lwork;
  xinfo = *info;

  dsygv_( &xitype, jobz, uplo, &xn,
          a, &xlda, b, &xldb,
          w, work, &xlwork, &xinfo );

 *info = (BlopexInt)xinfo;
}

BlopexInt zpotrf_redirect (char *uplo, BlopexInt *n, komplex *R, BlopexInt *lda, BlopexInt *info)
{
#ifdef BLOPEX_MATLAB_32
  int xn, xlda, xinfo;
#else
  mwSignedIndex xn, xlda, xinfo;
#endif
  xn = *n;
  xlda = *lda;
#ifdef BLOPEX_MATLAB_32
  zpotrf_(uplo, &xn, R, &xlda, &xinfo);
#else
  zpotrf_(uplo, &xn,(double *)R, &xlda, &xinfo);
#endif
  *info = (BlopexInt)xinfo;
}

BlopexInt zhegv_redirect (BlopexInt *itype, char *jobz, char *uplo, BlopexInt *
                    n, komplex *a, BlopexInt *lda, komplex *b, BlopexInt *ldb,
                    double *w, komplex *work, BlopexInt *lwork, double * rwork, BlopexInt *info)
{
#ifdef BLOPEX_MATLAB_32
  int xitype, xn, xlda, xldb, xlwork, xinfo;
#else
  mwSignedIndex xitype, xn, xlda, xldb, xlwork, xinfo;
#endif
  xitype = *itype;
  xn = *n;
  xlda = *lda;
  xldb = *ldb;
  xlwork = *lwork;
  xinfo = *info;
#ifdef BLOPEX_MATLAB_32
  zhegv_( &xitype, jobz, uplo, &xn,
          a, &xlda, b, &xldb,
          w, work, &xlwork, rwork, &xinfo );
#else
  zhegv_( &xitype, jobz, uplo, &xn,
          (double *)a, &xlda,(double *) b, &xldb,
          w, (double *)work, &xlwork, rwork, &xinfo );
#endif
 *info = (BlopexInt)xinfo;
}

int MATLABSetupInterpreter( mv_InterfaceInterpreter *i )
{
    /* Vector part */
  i->CreateVector = NULL;
  i->DestroyVector = NULL;
  i->InnerProd = NULL; 
  i->CopyVector = NULL;
  i->ClearVector = NULL;
  i->SetRandomValues = NULL;
  i->ScaleVector = NULL;
  i->Axpy = NULL;

  /* Multivector part */

  i->CreateMultiVector = NULL;
  i->CopyCreateMultiVector = CreateCopyMultiVector;
  i->DestroyMultiVector = DestroyMultiVector;

  i->Width = MultiVectorWidth;
  i->Height = NULL;
  i->SetMask = MultiSetMask;
  i->CopyMultiVector = CopyMultiVector;
  i->ClearMultiVector = NULL;
  i->SetRandomVectors = NULL;
  i->MultiInnerProd = MultiInnerProd;
  i->MultiInnerProdDiag = MultiInnerProdDiag;
  i->MultiVecMat = MultiVectorByMatrix;
  i->MultiVecMatDiag = MultiVectorByDiagonal;
  i->MultiAxpy = MultiVectorAxpy;
  i->MultiXapy = NULL;
  i->Eval = NULL;
 return 0;
}
