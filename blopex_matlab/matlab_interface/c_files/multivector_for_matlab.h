/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1) LGPL Version 2.1 or above.See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team http://code.google.com/p/blopex/     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* multivector structure for use with matlab */

#ifndef MULTIVECTOR_FOR_MATLAB
#define MULTIVECTOR_FOR_MATLAB

#include "matrix.h"
#include "fortran_matrix.h" /*eugene, add definition of complex data type*/

typedef struct matlabMultiVector* matlabMultiVectorPtr;

/* function prototypes follow */
matlabMultiVectorPtr MatlabMultiVectorCreate(const mxArray* blockVector);
void MatlabMultiVectorDestroy(matlabMultiVectorPtr X);
mxArray* MatlabMultiVectorGetCopyOfBlockVector(matlabMultiVectorPtr X);
BlopexInt MatlabMultiVectorWidth(matlabMultiVectorPtr X);
void MatlabMultiVectorSetMask(matlabMultiVectorPtr X, const BlopexInt* new_mask);
matlabMultiVectorPtr MatlabMultiVectorCopyCreate(matlabMultiVectorPtr X,
                                                 BlopexInt CopyValues);
void MatlabMultiVectorCopy(matlabMultiVectorPtr X, matlabMultiVectorPtr Y);
void MatlabMultiVectorInnerProd(matlabMultiVectorPtr X, 
                                matlabMultiVectorPtr Y, BlopexInt xyGHeight, 
                                BlopexInt xyHeight, BlopexInt xyWidth, void* xyVal); 
void MatlabMultiVectorInnerProdDiag(matlabMultiVectorPtr X,
                               matlabMultiVectorPtr Y, BlopexInt* mask, BlopexInt n, 
                               void* diag);                           
void MatlabMultiVectorByMatrix(matlabMultiVectorPtr X, BlopexInt rGHeight, 
                               BlopexInt rHeight, BlopexInt rWidth, void* rVal,  
                               matlabMultiVectorPtr Y);
void MatlabMultiVectorByDiag(matlabMultiVectorPtr X, BlopexInt *mask, BlopexInt n,
                               void *diag, matlabMultiVectorPtr Y);   
void MatlabMultiVectorAxpy(double alpha, matlabMultiVectorPtr X, 
                           matlabMultiVectorPtr Y);
void MatlabMultiVectorMatMultiVec( mxArray* data, matlabMultiVectorPtr X,
                                   matlabMultiVectorPtr Y );
#endif /* #ifndef MULTIVECTOR_FOR_MATLAB */
