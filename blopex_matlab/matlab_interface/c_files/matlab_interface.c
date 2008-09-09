/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* this file contains interface routines to matlab for blopex */
#include "matlab_interface.h"

#include "interpreter.h"

#include "multivector_for_matlab.h"

void*
CreateCopyMultiVector( void* src_, int copyValues )
{
  return MatlabMultiVectorCopyCreate((matlabMultiVectorPtr) src_,
                                     copyValues);
}

void
DestroyMultiVector( void *vvector )
{
  MatlabMultiVectorDestroy((matlabMultiVectorPtr) vvector);
}


int
MultiVectorWidth( void* v )
{
  return MatlabMultiVectorWidth((matlabMultiVectorPtr) v);
}

void
MultiSetMask( void *vector, int *mask )
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
				    int gh, int h, int w, double* v )
{
  MatlabMultiVectorInnerProd((matlabMultiVectorPtr) x_, 
                             (matlabMultiVectorPtr) y_, gh, h, w, v);
}

void
MultiInnerProdDiag( void* x_, void* y_,
					int* mask, int n, double* diag )
{
  MatlabMultiVectorInnerProdDiag((matlabMultiVectorPtr) x_,
                                 (matlabMultiVectorPtr) y_, mask, n, diag);
}


void
MultiVectorByDiagonal( void* x, 
				      int* mask, int n, double* diag,
				      void* y )
{
  MatlabMultiVectorByDiag((matlabMultiVectorPtr) x, mask, n, diag, 
                          (matlabMultiVectorPtr) y);
}


void 
MultiVectorByMatrix( void* x, 
			       int gh, int h, int w, double* v,
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
