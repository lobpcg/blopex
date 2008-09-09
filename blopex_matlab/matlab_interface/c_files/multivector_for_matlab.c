/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* multivector for matlab: structure and access routines */

#include "multivector_for_matlab.h"

#include "matrix.h"
#include "mex.h"

#include <stdlib.h>

struct matlabMultiVector
{
  mxArray* blockVector;        /* the pointer to the "blockVector" */
  mxArray* mask;               /* the pointer to the mask */
} ;

matlabMultiVectorPtr MatlabMultiVectorCreate(const mxArray* blockVector)
{
/* creates multivector from given "blockVector"; this also creates a
 mask where all vectors are "active"; NOTE that a copy of blockVector is
 made during creation of multivector */
  matlabMultiVectorPtr  X;
  mxArray*              mask;
  const int*            dimensions;
  mxLogical*            mask_data;
  int                   i;
  
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

int MatlabMultiVectorWidth(matlabMultiVectorPtr X)
{
  const int * dimensions = mxGetDimensions(X->blockVector);
  
  return(dimensions[1]);
}

void MatlabMultiVectorSetMask(matlabMultiVectorPtr X, const int* new_mask)
{
  mxLogical*  mask_data = mxGetLogicals(X->mask);
  int         i;
  int         mask_size = mxGetM(X->mask);
  
  if (new_mask==NULL)
    for (i=0; i<mask_size; i++)
      mask_data[i] = 1;
  else
    for (i=0; i<mask_size; i++)
      mask_data[i] = new_mask[i]? 1 : 0;
}

matlabMultiVectorPtr MatlabMultiVectorCopyCreate(matlabMultiVectorPtr X,
                                                 int CopyValues)
/* makes copy of multivector X ("blockVector" is copied, mask is set to all
   "active"); Parameter "CopyValues" is ignored for now */
{
  matlabMultiVectorPtr  Y;
  mxArray*              blockVector;
  mxArray*              mask;
  const int*            dimensions;
  mxLogical*            mask_data;
  int                   i;
  

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
                                matlabMultiVectorPtr Y, int xyGHeight, 
                                int xyHeight, int xyWidth, double* xyVal)
/* X(:,maskX)'*Y(:,maskY) is calculated and stored in fortran-style 
   matrix */
{
  mxArray*     matrix = NULL;
  mxArray*     prhs[] = {X->blockVector, X->mask, Y->blockVector, Y->mask};
  double*      matrix_data;
  int          gap = xyGHeight - xyHeight;
  int          i, j;
  
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

void MatlabMultiVectorInnerProdDiag(matlabMultiVectorPtr X,
                               matlabMultiVectorPtr Y, int* mask, int n, 
                               double* diag)
/* calculates diag(X(:,maskX)'*Y(:,maskY)) and stores it in unmasked 
   entries of "diag". X and Y mast have the same number of active vectors. 
   The same number of entries in diag should be unmasked */
{
  mxArray*     products_vector = NULL;
  double*      products_data;
  mxArray*     prhs[] = {X->blockVector, X->mask, Y->blockVector, Y->mask};
  int*         diag_active_ind;
  int          diag_num_active;
  int          i;
  
  mexCallMATLAB(1, &products_vector, 4, prhs,"matlabMultiInnerProdDiag");
  
  /* build list of active indices in diag */ 
  diag_active_ind = mxCalloc(n, sizeof(int));
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
  for (i=0; i<diag_num_active; i++)
    diag[diag_active_ind[i]] = products_data[i];
    
  mxFree(diag_active_ind);
  mxDestroyArray(products_vector);
}

void MatlabMultiVectorByMatrix(matlabMultiVectorPtr X, int rGHeight, 
                               int rHeight, int rWidth, double* rVal, 
                               matlabMultiVectorPtr Y)
/* Y(:,maskY)=X(:,maskX)*r; "r" is stored in fortran style */
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        mxCreateDoubleMatrix(rHeight, rWidth, mxREAL)};
  int        i,j;
  int        gap;
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

void MatlabMultiVectorByDiag(matlabMultiVectorPtr X, int *mask, int n,
                             double *diag, matlabMultiVectorPtr Y)
/* Y(:,maskY)=X(:,maskX)*diag(alpha(mask)) */
{
  mxArray*   new_blockVectorY;
  mxArray*   prhs[5] = {X->blockVector, X->mask, Y->blockVector, Y->mask,
                        NULL};
  int*       diag_active_ind;
  int        diag_num_active;
  int        i;
  double*    matrix_data;
  
    /* build list of active indices in diag */ 
  diag_active_ind = mxCalloc(n, sizeof(int));
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
