/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ University of Colorado at Denver                                  */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* multivector structure for use with matlab */
#ifndef MULTIVECTOR_FOR_MATLAB
#define MULTIVECTOR_FOR_MATLAB

#include "matrix.h"

typedef struct matlabMultiVector* matlabMultiVectorPtr;

/* function prototypes follow */
matlabMultiVectorPtr MatlabMultiVectorCreate(const mxArray* blockVector);
void MatlabMultiVectorDestroy(matlabMultiVectorPtr X);
mxArray* MatlabMultiVectorGetCopyOfBlockVector(matlabMultiVectorPtr X);
int MatlabMultiVectorWidth(matlabMultiVectorPtr X);
void MatlabMultiVectorSetMask(matlabMultiVectorPtr X, const int* new_mask);
matlabMultiVectorPtr MatlabMultiVectorCopyCreate(matlabMultiVectorPtr X,
                                                 int CopyValues);
void MatlabMultiVectorCopy(matlabMultiVectorPtr X, matlabMultiVectorPtr Y);
void MatlabMultiVectorInnerProd(matlabMultiVectorPtr X, 
                                matlabMultiVectorPtr Y, int xyGHeight, 
                                int xyHeight, int xyWidth, double* xyVal);
void MatlabMultiVectorInnerProdDiag(matlabMultiVectorPtr X,
                               matlabMultiVectorPtr Y, int* mask, int n, 
                               double* diag);
void MatlabMultiVectorByMatrix(matlabMultiVectorPtr X, int rGHeight, 
                               int rHeight, int rWidth, double* rVal, 
                               matlabMultiVectorPtr Y);
void MatlabMultiVectorByDiag(matlabMultiVectorPtr X, int *mask, int n,
                             double *diag, matlabMultiVectorPtr Y);
void MatlabMultiVectorAxpy(double alpha, matlabMultiVectorPtr X, 
                           matlabMultiVectorPtr Y);
void MatlabMultiVectorMatMultiVec( mxArray* data, matlabMultiVectorPtr X,
                                   matlabMultiVectorPtr Y );
#endif /* #ifndef MULTIVECTOR_FOR_MATLAB */
