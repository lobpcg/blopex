/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.2) MIT / Apache-2.0 License                   */
/* @@@ Copyright 2010-2019 Team https://github.com/lobpcg/blopex       */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

typedef struct {double real, imag;} komplex;

#include "multi_vector.h"
#include "matmultivec.h"
#include "pcg_multi.h"
#include "interpreter.h"
#include <assert.h>
#include <string.h>

/*--------------------------------------------------------------------------
  CreateCopyMultiVector                                              generic
  --------------------------------------------------------------------------*/
void*
CreateCopyMultiVector( void* src_, BlopexInt copyValues )
{
   serial_Multi_Vector *src = (serial_Multi_Vector *)src_;
   serial_Multi_Vector *dest;

   /* create vector with the same parameters as src */

   dest = serial_Multi_VectorCreate(src->size, src->num_vectors);
   serial_Multi_VectorInitialize(dest);

   /* copy values if necessary */

   if (copyValues)
      serial_Multi_VectorCopyWithoutMask(src, dest);

   return dest;
}

/*--------------------------------------------------------------------------
  DestroyMultiVector                                                 generic
  --------------------------------------------------------------------------*/
void
DestroyMultiVector( void *vvector )
{
  BlopexInt dummy;
  serial_Multi_Vector *vector = vvector;

  dummy=serial_Multi_VectorDestroy( vector );
}

/*--------------------------------------------------------------------------
  MultiVectorWidth                                                   generic
  --------------------------------------------------------------------------*/
BlopexInt
MultiVectorWidth( void* v )
{
  return ((serial_Multi_Vector*)v)->num_vectors;
}

/*--------------------------------------------------------------------------
  MultiSetMask                                                       generic
  --------------------------------------------------------------------------*/
void
MultiSetMask( void *vector, BlopexInt *mask )
{
   serial_Multi_VectorSetMask( ( serial_Multi_Vector *)vector, mask );
}

/*--------------------------------------------------------------------------
  CopyMultiVector                                                    complex
  --------------------------------------------------------------------------*/
void
CopyMultiVector( void *x, void *y)
{
   BlopexInt dummy;

   dummy = serial_Multi_VectorCopy( (serial_Multi_Vector *) x,
                                    (serial_Multi_Vector *) y);
}

/*--------------------------------------------------------------------------
  ClearMultiVector                                                   complex
  --------------------------------------------------------------------------*/
void
ClearMultiVector(void *x)
{
   komplex kzero = {0.0,0.0};
   BlopexInt dummy;

   dummy=serial_Multi_VectorSetConstantValues((serial_Multi_Vector *)x,kzero);
}

/*--------------------------------------------------------------------------
  MultiVectorSetRandomValues                                         complex
  --------------------------------------------------------------------------*/
void
SetMultiVectorRandomValues(void *x, BlopexInt seed)
{
   BlopexInt dummy;

   dummy= serial_Multi_VectorSetRandomValues((serial_Multi_Vector *) x, seed) ;
}


/*--------------------------------------------------------------------------
  MultiInnerProd                                                      complex
  --------------------------------------------------------------------------*/
void
MultiInnerProd(void * x_, void * y_,
                    BlopexInt gh, BlopexInt h, BlopexInt w, void* v )
{
   serial_Multi_VectorInnerProd( (serial_Multi_Vector *)x_,
                                 (serial_Multi_Vector *)y_,
                                 gh, h, w, (komplex *) v);
}


/*--------------------------------------------------------------------------
  MultiInnerProdDiag                                                  complex
  --------------------------------------------------------------------------*/
void
MultiInnerProdDiag( void* x_, void* y_,
                    BlopexInt* mask, BlopexInt n, void* diag )
{
   serial_Multi_VectorInnerProdDiag( (serial_Multi_Vector *)x_,
                                     (serial_Multi_Vector *)y_,
                                      mask, n, (komplex *) diag);
}

/*--------------------------------------------------------------------------
  MultiVectorByDiagonal                                              complex
  --------------------------------------------------------------------------*/
void
MultiVectorByDiagonal( void* x,
                       BlopexInt* mask, BlopexInt n, void* diag,
                       void* y )
{
   BlopexInt dummy;

   dummy = serial_Multi_VectorByDiag( (serial_Multi_Vector *) x, mask, n,
                                      (komplex *) diag,
                                      (serial_Multi_Vector *) y );
}

/*--------------------------------------------------------------------------
  MultiVectorByMatrix                                                complex
  --------------------------------------------------------------------------*/
void
MultiVectorByMatrix( void* x,
                   BlopexInt gh, BlopexInt h, BlopexInt w, void* v,
                   void* y )
{
   serial_Multi_VectorByMatrix((serial_Multi_Vector *)x, gh, h,
                                w, (komplex *) v, (serial_Multi_Vector *)y);

}
/*--------------------------------------------------------------------------
  MultiAxpy                                                          complex
  --------------------------------------------------------------------------*/
void
MultiVectorAxpy( double alpha, void   *x, void   *y)
{
   serial_Multi_VectorAxpy(  alpha,
                            (serial_Multi_Vector *) x,
                            (serial_Multi_Vector *) y) ;
}
/*--------------------------------------------------------------------------
  MultiVectorPrint                                                   complex
  --------------------------------------------------------------------------*/
void
MultiVectorPrint( void *x, char * tag, BlopexInt limit )
{
   serial_Multi_VectorPrint( (serial_Multi_Vector *) x, tag, limit );
}
/*--------------------------------------------------------------------------
  MatMultiVec                                                        complex
  --------------------------------------------------------------------------*/
void MatMultiVec (void * A, void * X, void * AX)
{
   /* A,X,AX are pointers to serial_Multi_Vector's
      if you are here then blopex_solve was called with operatorA = MatMultiVec

      A may be NULL if blopex_solve called with operatorA only
      and no matrix operatorAData, in that case a call to the actual operator
      must be done, i.e. 1D Laplacian

      If A is not NULL then compute AX = A*X */

   serial_Multi_Vector * sA;
   serial_Multi_Vector * sX;
   serial_Multi_Vector * sAX;

   /* 1D Laplacian operator, Dirichlet BC */

   if (A==NULL)
     serial_Lapl1DMatMultiVec( (serial_Multi_Vector *) X,
                               (serial_Multi_Vector *) AX  );
   else {
     sA  = (serial_Multi_Vector *)A;
     sX  = (serial_Multi_Vector *)X;
     sAX = (serial_Multi_Vector *)AX;
     serial_Multi_VectorByMulti_Vector(sA, sX, sAX);
   }
}

BlopexInt
/*--------------------------------------------------------------------------
  SerialSetupInterpreter                                             complex
  --------------------------------------------------------------------------*/
SerialSetupInterpreter( mv_InterfaceInterpreter *i )
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
  i->ClearMultiVector = ClearMultiVector;
  i->SetRandomVectors = SetMultiVectorRandomValues;
  i->MultiInnerProd = MultiInnerProd;
  i->MultiInnerProdDiag = MultiInnerProdDiag;
  i->MultiVecMat = MultiVectorByMatrix;
  i->MultiVecMatDiag = MultiVectorByDiagonal;
  i->MultiAxpy = MultiVectorAxpy;
  i->MultiXapy = NULL;
  i->Eval = NULL;
  i->MultiPrint = MultiVectorPrint;

  return 0;
}
