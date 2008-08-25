/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#include "multi_vector.h"
#include "matmultivec.h"
#include "pcg_multi.h"
#include "interpreter.h"
#include <assert.h>
#include <string.h>

void*
CreateCopyMultiVector( void* src_, int copyValues )
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
 * DestroyMultiVector
 *--------------------------------------------------------------------------*/

void
DestroyMultiVector( void *vvector )
{
  int dummy; 
  serial_Multi_Vector *vector = vvector;

  dummy=serial_Multi_VectorDestroy( vector );
}


int
MultiVectorWidth( void* v )
{
  return ((serial_Multi_Vector*)v)->num_vectors;
}


/*--------------------------------------------------------------------------
 * MultiSetMask
 *--------------------------------------------------------------------------*/

void
MultiSetMask( void *vector, int *mask )
{
   serial_Multi_VectorSetMask( ( serial_Multi_Vector *)vector, mask );
}


/*--------------------------------------------------------------------------
 * CopyMultiVector
 *--------------------------------------------------------------------------*/

void
CopyMultiVector( void *x, void *y)
{
   int dummy;
    
   dummy = serial_Multi_VectorCopy( (serial_Multi_Vector *) x,
                                      (serial_Multi_Vector *) y);
}

/*--------------------------------------------------------------------------
 * hypre_ParKrylovClearMultiVector
 *--------------------------------------------------------------------------*/

void
ClearMultiVector(void *x)
{
   int dummy;
   
   dummy=serial_Multi_VectorSetConstantValues((serial_Multi_Vector *)x,0.0);
}


/*--------------------------------------------------------------------------
 * MultiVectorSetRandomValues
 *--------------------------------------------------------------------------*/

void
SetMultiVectorRandomValues(void *x, int seed)
{
   int dummy;
   
   dummy= serial_Multi_VectorSetRandomValues((serial_Multi_Vector *) x, seed) ;
}


/*--------------------------------------------------------------------------
 * MultiInnerProd
 *--------------------------------------------------------------------------*/

void
MultiInnerProd(void * x_, void * y_,
				    int gh, int h, int w, double* v )
{
   serial_Multi_VectorInnerProd( (serial_Multi_Vector *)x_,
                                 (serial_Multi_Vector *)y_, 
                                 gh, h, w, v);
}



/*--------------------------------------------------------------------------
 * MultiInnerProdDiag
 *--------------------------------------------------------------------------*/

void
MultiInnerProdDiag( void* x_, void* y_,
					int* mask, int n, double* diag )
{
   serial_Multi_VectorInnerProdDiag( (serial_Multi_Vector *)x_,
                                     (serial_Multi_Vector *)y_, 
				     mask, n, diag);
}


void
MultiVectorByDiagonal( void* x, 
				      int* mask, int n, double* diag,
				      void* y )
{
   int dummy;
   
   dummy = serial_Multi_VectorByDiag( (serial_Multi_Vector *) x, mask, n, diag, 
                                                      (serial_Multi_Vector *) y );
}


void 
MultiVectorByMatrix( void* x, 
			       int gh, int h, int w, double* v,
			       void* y )
{                               
   serial_Multi_VectorByMatrix((serial_Multi_Vector *)x, gh, h, 
                                w, v, (serial_Multi_Vector *)y);

}
/*--------------------------------------------------------------------------
 * MultiAxpy
 *--------------------------------------------------------------------------*/

void
MultiVectorAxpy( double alpha, void   *x, void   *y)
{
   serial_Multi_VectorAxpy(  alpha, 
                              (serial_Multi_Vector *) x,
                              (serial_Multi_Vector *) y) ;
}

void MatMultiVec (void * data, void * x, void * y)
{
   /* temporary implementation: diagonal matrix */
   serial_Lapl1DMatMultiVec( (serial_Multi_Vector *) x,
                           (serial_Multi_Vector *) y  );
}

int
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

  return 0;
}
