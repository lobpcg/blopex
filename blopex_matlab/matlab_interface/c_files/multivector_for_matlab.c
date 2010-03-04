/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1) LGPL Version 3 or above.  See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team http://code.google.com/p/blopex/     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* multivector for matlab: structure and access routines */

extern int BLOPEX_MATLAB_COMPLEX;

#include "multivector_for_matlab.h"

#include "matrix.h"
#include "mex.h"
#include "fortran_matrix.h"                                                               

#include <stdlib.h>

struct matlabMultiVector
{
  mxArray* blockVector;        /* the pointer to the "blockVector" */
  mxArray* mask;               /* the pointer to the mask */
} ;

/* declare complex/real function prototypes */
/* ---- complex arithmetic ----*/
  void MatlabMultiVectorInnerProd_complex(matlabMultiVectorPtr,
                                          matlabMultiVectorPtr, BlopexInt,
                                          BlopexInt,BlopexInt, komplex*); 
  void MatlabMultiVectorInnerProdDiag_complex(matlabMultiVectorPtr,
                               matlabMultiVectorPtr, BlopexInt*, BlopexInt, 
                               komplex*);
  void MatlabMultiVectorByMatrix_complex(matlabMultiVectorPtr, BlopexInt, 
                               BlopexInt, BlopexInt, komplex*, 
                               matlabMultiVectorPtr);
  void MatlabMultiVectorByDiag_complex(matlabMultiVectorPtr, BlopexInt*, BlopexInt,
                               komplex*, matlabMultiVectorPtr);
 
/* ---- real arithmetic ---- */
  void MatlabMultiVectorInnerProd_double(matlabMultiVectorPtr,
                                matlabMultiVectorPtr, BlopexInt,
                                BlopexInt,BlopexInt, double*);
  void MatlabMultiVectorInnerProdDiag_double(matlabMultiVectorPtr,
                                matlabMultiVectorPtr, BlopexInt*, BlopexInt, 
                                double*);
  void MatlabMultiVectorByMatrix_double(matlabMultiVectorPtr, BlopexInt, 
                                BlopexInt, BlopexInt, double*, 
                                matlabMultiVectorPtr);
  void MatlabMultiVectorByDiag_double(matlabMultiVectorPtr, BlopexInt*, BlopexInt,
                                double*, matlabMultiVectorPtr);

matlabMultiVectorPtr MatlabMultiVectorCreate(const mxArray* blockVector)
{
/* creates multivector from given "blockVector"; this also creates a
 mask where all vectors are "active"; NOTE that a copy of blockVector is
 made during creation of multivector */
  matlabMultiVectorPtr  X;
  mxArray*              mask;
  const mwSize*         dimensions;
  mxLogical*            mask_data;
  BlopexInt                   i;
  
   /* get size of blockVector to build adequate mask */
  if (mxGetNumberOfDimensions(blockVector) !=2)
    mexErrMsgTxt("blockVector is NOT 2d array");
  dimensions = mxGetDimensions(blockVector);
  if (dimensions[0]<dimensions[1])
    mexErrMsgTxt("blockVector must be tall, not fat");
  mask = mxCreateLogicalMatrix(dimensions[1], 1);
  mask_data = mxGetLogicals(mask);
  for (i=0; i<dimensions[1];i++)
    mask_data[i] = 1;
  
  X = mxCalloc(1, sizeof(struct matlabMultiVector));
  X->blockVector = mxDuplicateArray(blockVector);
  X->mask = mask;
  return X;
}

void MatlabMultiVectorDestroy(matlabMultiVectorPtr X)
{
  mxDestroyArray(X->blockVector);
  mxDestroyArray(X->mask);
  mxFree(X);
}

mxArray* MatlabMultiVectorGetCopyOfBlockVector(matlabMultiVectorPtr X)
{
  return (mxDuplicateArray(X->blockVector));
}

BlopexInt MatlabMultiVectorWidth(matlabMultiVectorPtr X)
{
  const mwSize * dimensions = mxGetDimensions(X->blockVector);
  return(dimensions[1]);
}

void MatlabMultiVectorSetMask(matlabMultiVectorPtr X, const BlopexInt* new_mask)
{
  mxLogical*  mask_data = mxGetLogicals(X->mask);
  BlopexInt         i;
  BlopexInt         mask_size = mxGetM(X->mask);
  
  if (new_mask==NULL)
    for (i=0; i<mask_size; i++)
      mask_data[i] = 1;
  else
    for (i=0; i<mask_size; i++)
      mask_data[i] = new_mask[i]? 1 : 0;
}

matlabMultiVectorPtr MatlabMultiVectorCopyCreate(matlabMultiVectorPtr X,
                                                 BlopexInt CopyValues)
/* makes copy of multivector X ("blockVector" is copied, mask is set to all
   "active"); Parameter "CopyValues" is ignored for now */
{
  matlabMultiVectorPtr  Y;
  mxArray*              blockVector;
  mxArray*              mask;
  const mwSize*            dimensions;
  mxLogical*            mask_data;
  BlopexInt                   i;
  

  blockVector = mxDuplicateArray(X->blockVector);
  
  /* get size of blockVector to build adequate mask */
  dimensions = mxGetDimensions(blockVector);
  mask = mxCreateLogicalMatrix(dimensions[1], 1);
  mask_data = mxGetLogicals(mask);
  for (i=0; i<dimensions[1];i++)
    mask_data[i] = 1;
  
  Y = mxCalloc(1, sizeof(struct matlabMultiVector));
  Y->blockVector = blockVector;
  Y->mask = mask;
  return Y;
}

void MatlabMultiVectorCopy(matlabMultiVectorPtr X, matlabMultiVectorPtr Y)
/* copy multivector "x" into "y" with respect to a mask */
{
  mxArray*     lhs = NULL;
  mxArray*     prhs[] = {X->blockVector, X->mask, Y->blockVector, Y->mask};
  
  mexCallMATLAB(1, &lhs, 4, prhs, "matlabCopyMultiVector");
  mxDestroyArray(Y->blockVector);
  Y->blockVector = lhs;
}

void MatlabMultiVectorInnerProd(matlabMultiVectorPtr X, 
                                matlabMultiVectorPtr Y, BlopexInt xyGHeight, 
                                BlopexInt xyHeight, BlopexInt xyWidth, void* xyVal)
/* X(:,maskX)'*Y(:,maskY) is calculated and stored in fortran-style 
   matrix */
{
  if (BLOPEX_MATLAB_COMPLEX) 
    MatlabMultiVectorInnerProd_complex(X,Y,xyGHeight,xyHeight,xyWidth,(komplex*)xyVal);
  else
    MatlabMultiVectorInnerProd_double(X,Y,xyGHeight,xyHeight,xyWidth,(double*)xyVal);
}

void MatlabMultiVectorInnerProdDiag(matlabMultiVectorPtr X,
                               matlabMultiVectorPtr Y, BlopexInt* mask, BlopexInt n, 
                               void* diag)
/* calculates diag(X(:,maskX)'*Y(:,maskY)) and stores it in unmasked 
   entries of "diag". X and Y mast have the same number of active vectors. 
   The same number of entries in diag should be unmasked */
{
  if (BLOPEX_MATLAB_COMPLEX)
    MatlabMultiVectorInnerProdDiag_complex(X,Y,mask,n,(komplex*)diag);
  else
    MatlabMultiVectorInnerProdDiag_double(X,Y,mask,n,(double*)diag);
}

void MatlabMultiVectorByMatrix(matlabMultiVectorPtr X, BlopexInt rGHeight, 
                               BlopexInt rHeight, BlopexInt rWidth, void* rVal, 
                               matlabMultiVectorPtr Y)
/* Y(:,maskY)=X(:,maskX)*r; "r" is stored in fortran style */
{
  if (BLOPEX_MATLAB_COMPLEX)
    MatlabMultiVectorByMatrix_complex(X,rGHeight,rHeight,rWidth,(komplex*)rVal,Y);
  else
    MatlabMultiVectorByMatrix_double(X,rGHeight,rHeight,rWidth,(double*)rVal,Y);
}

void MatlabMultiVectorByDiag(matlabMultiVectorPtr X, BlopexInt *mask, BlopexInt n,
                             void *diag, matlabMultiVectorPtr Y)
/* Y(:,maskY)=X(:,maskX)*diag(alpha(mask)) */
{
  if (BLOPEX_MATLAB_COMPLEX)
    MatlabMultiVectorByDiag_complex(X,mask,n,(komplex*)diag,Y);
  else
    MatlabMultiVectorByDiag_double(X,mask,n,(double*)diag,Y);
}

void MatlabMultiVectorAxpy(double alpha, matlabMultiVectorPtr X,
                          matlabMultiVectorPtr Y)
/* Y(:,maskY) = alpha*X(:,maskX)+Y(:,maskY) */
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        mxCreateDoubleScalar(alpha)};
                        
  mexCallMATLAB(1, &new_blockVectorY, 5, prhs,"matlabMultiAxpy");
  
  mxDestroyArray(Y->blockVector);
  Y->blockVector = new_blockVectorY;

  mxDestroyArray(prhs[4]);
}

void MatlabMultiVectorMatMultiVec( mxArray* data, matlabMultiVectorPtr X,
                                   matlabMultiVectorPtr Y )
/* Y(:,maskY) = operator * X(:,maskX) */
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        data};
                        
  mexCallMATLAB(1, &new_blockVectorY, 5, prhs,"matlabMultiMatVec");
  
  mxDestroyArray(Y->blockVector);
  Y->blockVector = new_blockVectorY;
}

/* --------- implementation of complex/real matlab multivector operations -------- */
/* ---- complex arithmetic ----*/    
void MatlabMultiVectorInnerProd_complex(matlabMultiVectorPtr X,
                                        matlabMultiVectorPtr Y, BlopexInt xyGHeight,
                                        BlopexInt xyHeight,BlopexInt xyWidth, komplex* xyVal)
{
  mxArray*     matrix = NULL;
  mxArray*     prhs[] = {X->blockVector, X->mask, Y->blockVector, Y->mask};
  double*      matrix_data_real;
  double*      matrix_data_imag;
  komplex      temp = {0,0};
  BlopexInt          gap = xyGHeight - xyHeight;
  BlopexInt          i, j;
  
  mexCallMATLAB(1, &matrix, 4, prhs,"matlabMultiInnerProd");
  
  matrix_data_real = mxGetPr(matrix);
  matrix_data_imag = mxGetPi(matrix);
    
  mxAssert(matrix_data_real != NULL || matrix_data_imag != NULL,
            "matrix data is NULL pointer");
  
  mxAssert(mxGetM(matrix)==xyHeight && mxGetN(matrix)==xyWidth,
           "inconsistent matrix-vector geometry");
  
  for(j=0; j<xyWidth; j++)
  {
    for(i=0; i<xyHeight; i++)
    {
      if (matrix_data_real != NULL) temp.real = *(matrix_data_real++);
      if (matrix_data_imag != NULL) temp.imag = *(matrix_data_imag++);
      *xyVal++ = temp;
      temp.real = 0;
      temp.imag = 0;
    }
    xyVal += gap;
  }
  
  mxDestroyArray(matrix);
}

void MatlabMultiVectorInnerProdDiag_complex(matlabMultiVectorPtr X,
                               matlabMultiVectorPtr Y, BlopexInt* mask, BlopexInt n, 
                               komplex* diag)
{
  mxArray*     products_vector = NULL;
  double*      products_data_real;
  double*      products_data_imag;
  mxArray*     prhs[] = {X->blockVector, X->mask, Y->blockVector, Y->mask};
  BlopexInt*         diag_active_ind;
  BlopexInt          diag_num_active;
  BlopexInt          i;
  komplex      temp = {0,0};
  
  mexCallMATLAB(1, &products_vector, 4, prhs,"matlabMultiInnerProdDiag");
  
  /* build list of active indices in diag */ 
  diag_active_ind = mxCalloc(n, sizeof(BlopexInt));
  diag_num_active = 0;
  
  if (mask!=NULL)
    for (i=0; i<n; i++)
  {
    if (mask[i])
      diag_active_ind[diag_num_active++]=i;
    }
  else
    for (i=0; i<n; i++)
      diag_active_ind[diag_num_active++]=i;
  
  mxAssert(diag_num_active==mxGetM(products_vector) && 
           mxGetN(products_vector)==1,"wrong matrix-vector-mask geometry");
  
  /* copy calculated inner products to unmasked entries of "diag" */
  products_data_real = mxGetPr(products_vector);
  products_data_imag = mxGetPi(products_vector);
  mxAssert((products_data_real!=NULL) || (products_data_imag != NULL),"products_data is NULL pointer");
  for (i=0; i<diag_num_active; i++)
  {
    if (products_data_real != NULL) temp.real = products_data_real[i];
    if (products_data_imag != NULL) temp.imag = products_data_imag[i];
    diag[diag_active_ind[i]] = temp;
    temp.real = 0;
    temp.imag = 0;
  }  
  mxFree(diag_active_ind);
  mxDestroyArray(products_vector);
}

void MatlabMultiVectorByMatrix_complex(matlabMultiVectorPtr X, BlopexInt rGHeight, 
                               BlopexInt rHeight, BlopexInt rWidth, komplex* rVal, 
                               matlabMultiVectorPtr Y)
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        mxCreateDoubleMatrix(rHeight, rWidth, mxCOMPLEX)};
  BlopexInt        i,j;
  BlopexInt        gap;
  double*    matrix_data_real;
  double*    matrix_data_imag;
  komplex    temp;
 
  /* copy data from rVal to matlab matrix */
  matrix_data_real = mxGetPr(prhs[4]);
  matrix_data_imag = mxGetPi(prhs[4]);
  gap = rGHeight - rHeight;
  for(j=0; j<rWidth; j++)
  {
    for(i=0; i<rHeight; i++)
    {
      temp = *(rVal++);
      *(matrix_data_real++) = temp.real;
      *(matrix_data_imag++) = temp.imag;
    }
    rVal += gap;
  }
    
  mexCallMATLAB(1, &new_blockVectorY, 5, prhs,"matlabMultiVectorByMatrix");
  
  mxDestroyArray(Y->blockVector);
  Y->blockVector = new_blockVectorY;
  
  mxDestroyArray(prhs[4]);
}

void MatlabMultiVectorByDiag_complex(matlabMultiVectorPtr X, BlopexInt *mask, BlopexInt n,
                             komplex *diag, matlabMultiVectorPtr Y)
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        NULL};
  BlopexInt*       diag_active_ind;
  BlopexInt        diag_num_active;
  BlopexInt        i;
  double*    matrix_data_real;
  double*    matrix_data_imag;
  
    /* build list of active indices in diag */ 
  diag_active_ind = mxCalloc(n, sizeof(BlopexInt));
  diag_num_active = 0;
  if (mask!=NULL)
    for (i=0; i<n; i++)
  {
    if (mask[i])
      diag_active_ind[diag_num_active++]=i;
    }
  else
    for (i=0; i<n; i++)
      diag_active_ind[diag_num_active++]=i;
  
  /* copy values from active entries of diag to matlab array */
  prhs[4] = mxCreateDoubleMatrix(diag_num_active, 1, mxCOMPLEX);
  matrix_data_real = mxGetPr(prhs[4]);
  matrix_data_imag = mxGetPi(prhs[4]);

  for (i=0; i<diag_num_active; i++)
  {
    matrix_data_real[i] = diag[diag_active_ind[i]].real;
    matrix_data_imag[i] = diag[diag_active_ind[i]].imag;
  }
  mexCallMATLAB(1, &new_blockVectorY, 5, prhs,"matlabMultiVectorByDiag");
  
  mxDestroyArray(Y->blockVector);
  Y->blockVector = new_blockVectorY;
  
  mxFree(diag_active_ind);
  mxDestroyArray(prhs[4]);
}


/* ---- real arithmetic ---- */

void MatlabMultiVectorInnerProd_double(matlabMultiVectorPtr X,
                                       matlabMultiVectorPtr Y, BlopexInt xyGHeight,
                                       BlopexInt xyHeight,BlopexInt xyWidth, double* xyVal)
{
  mxArray*     matrix = NULL;
  mxArray*     prhs[] = {X->blockVector, X->mask, Y->blockVector, Y->mask};
  double*      matrix_data;
  BlopexInt          gap = xyGHeight - xyHeight;
  BlopexInt          i, j;
  
  mexCallMATLAB(1, &matrix, 4, prhs,"matlabMultiInnerProd");
  
  /* now copy data from "matrix" to "xyVal", then destroy "matrix" */
  if (mxIsComplex(matrix))
    mexErrMsgTxt("complex arithmetic not supported");
  
  matrix_data = mxGetPr(matrix);
  mxAssert((matrix_data!=NULL),"matrix data is NULL pointer");
  
  mxAssert(mxGetM(matrix)==xyHeight && mxGetN(matrix)==xyWidth,
           "inconsistent matrix-vector geometry");
  
  for(j=0; j<xyWidth; j++)
  {
    for(i=0; i<xyHeight; i++)
      *(xyVal++) = *(matrix_data++);
    xyVal += gap;
  }
  
  mxDestroyArray(matrix);
}

void MatlabMultiVectorInnerProdDiag_double(matlabMultiVectorPtr X,
                               matlabMultiVectorPtr Y, BlopexInt* mask, BlopexInt n, 
                               double* diag)
{
  mxArray*     products_vector = NULL;
  double*      products_data;
  mxArray*     prhs[] = {X->blockVector, X->mask, Y->blockVector, Y->mask};
  BlopexInt*         diag_active_ind;
  BlopexInt          diag_num_active;
  BlopexInt          i;
  
  mexCallMATLAB(1, &products_vector, 4, prhs,"matlabMultiInnerProdDiag");
  
  /* build list of active indices in diag */ 
  diag_active_ind = mxCalloc(n, sizeof(BlopexInt));
  diag_num_active = 0;
  
  if (mask!=NULL)
    for (i=0; i<n; i++)
  {
    if (mask[i])
      diag_active_ind[diag_num_active++]=i;
    }
  else
    for (i=0; i<n; i++)
      diag_active_ind[diag_num_active++]=i;
  
  mxAssert(diag_num_active==mxGetM(products_vector) && 
           mxGetN(products_vector)==1,"wrong matrix-vector-mask geometry");
  if (mxIsComplex(products_vector))
    mexErrMsgTxt("complex arithmetic not supported");  
  
  /* copy calculated inner products to unmasked entries of "diag" */
  products_data = mxGetPr(products_vector);
  mxAssert((products_data!=NULL),"products_data is NULL pointer");
  for (i=0; i<diag_num_active; i++)
    diag[diag_active_ind[i]] = products_data[i];
    
  mxFree(diag_active_ind);
  mxDestroyArray(products_vector);
}

void MatlabMultiVectorByMatrix_double(matlabMultiVectorPtr X, BlopexInt rGHeight, 
                               BlopexInt rHeight, BlopexInt rWidth, double* rVal, 
                               matlabMultiVectorPtr Y)
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        mxCreateDoubleMatrix(rHeight, rWidth, mxREAL)};
  BlopexInt        i,j;
  BlopexInt        gap;
  double*    matrix_data;
 
  /* copy data from rVal to matlab matrix */
  matrix_data = mxGetPr(prhs[4]);
  gap = rGHeight - rHeight;
  for(j=0; j<rWidth; j++)
  {
    for(i=0; i<rHeight; i++)
      *(matrix_data++) = *(rVal++);
    rVal += gap;
  }
    
  mexCallMATLAB(1, &new_blockVectorY, 5, prhs,"matlabMultiVectorByMatrix");
  
  mxDestroyArray(Y->blockVector);
  Y->blockVector = new_blockVectorY;
  
  mxDestroyArray(prhs[4]);
}

void MatlabMultiVectorByDiag_double(matlabMultiVectorPtr X, BlopexInt *mask, BlopexInt n,
                             double *diag, matlabMultiVectorPtr Y)
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        NULL};
  BlopexInt*       diag_active_ind;
  BlopexInt        diag_num_active;
  BlopexInt        i;
  double*    matrix_data;
  
    /* build list of active indices in diag */ 
  diag_active_ind = mxCalloc(n, sizeof(BlopexInt));
  diag_num_active = 0;
  if (mask!=NULL)
    for (i=0; i<n; i++)
  {
    if (mask[i])
      diag_active_ind[diag_num_active++]=i;
    }
  else
    for (i=0; i<n; i++)
      diag_active_ind[diag_num_active++]=i;
  
  /* copy values from active entries of diag to matlab array */
  prhs[4] = mxCreateDoubleMatrix(diag_num_active, 1, mxREAL);
  matrix_data = mxGetPr(prhs[4]);
  for (i=0; i<diag_num_active; i++)
    matrix_data[i] = diag[diag_active_ind[i]] ;
  
  mexCallMATLAB(1, &new_blockVectorY, 5, prhs,"matlabMultiVectorByDiag");
  
  mxDestroyArray(Y->blockVector);
  Y->blockVector = new_blockVectorY;
  
  mxFree(diag_active_ind);
  mxDestroyArray(prhs[4]);
}

