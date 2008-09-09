/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#include "mex.h"
#include "matrix.h"

#include "multivector_for_matlab.h"
#include "matlab_interface.h"

#include "lobpcg.h"

#define REQUIRED_NUM_INPUTS    8
#define REQUIRED_NUM_OUTPUTS   6

int dsygv_ (int *itype, char *jobz, char *uplo, int *
                    n, double *a, int *lda, double *b, int *ldb,
                    double *w, double *work, int *lwork, int *info);

int dpotrf_ (char *uplo, int *n, double *a, int *
                    lda, int *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  matlabMultiVectorPtr         x;
  mv_MultiVectorPtr            xx;
  matlabMultiVectorPtr         y;
  mv_MultiVectorPtr            yy;
  lobpcg_BLASLAPACKFunctions   blap_fn;
  int                          iterations;
  lobpcg_Tolerance             lobpcg_tol;
  mv_InterfaceInterpreter      ii;
  double*                      eigs;
  double*                      resid;
  double*                      eigs_hist;
  double*                      resid_hist;
  int                          maxit;
  int                          verbosity_level;
  int                          width;
  int                          exitflag;
  
/* Check for proper number of arguments. */
  if (nrhs != REQUIRED_NUM_INPUTS)
  {
    mexErrMsgTxt("wrong number of inputs");
  }
  else if (nlhs != REQUIRED_NUM_OUTPUTS)
  {
    mexErrMsgTxt("wrong number of outputs");
  }
  
  /* retrieve some input arguments */
  maxit = mxGetScalar(prhs[6]);
  lobpcg_tol.absolute = mxGetScalar(prhs[5]);
  lobpcg_tol.relative = 1e-50;
  verbosity_level = mxGetScalar(prhs[7]);
  width = mxGetN(prhs[0]);
  
  /* create multivector */
  x = MatlabMultiVectorCreate(prhs[0]);
  

 /* get memory for eigenvalues, eigenvalue history, residual norms,
     residual norms history */
  
   /* request memory for eig-vals */
  plhs[1] = mxCreateDoubleMatrix(width, 1, mxREAL);
  eigs = mxGetPr(plhs[1]);
  
   /* request mem. for eig-history */
  plhs[3] = mxCreateDoubleMatrix(width, maxit+1, mxREAL);
  eigs_hist = mxGetPr(plhs[3]);
  
   /* request memory for resid. norms */
  resid = mxCalloc(width, sizeof(double));
  
  /* req. mem. for resid-hist */
  plhs[4] = mxCreateDoubleMatrix(width, maxit+1, mxREAL);
  resid_hist = mxGetPr(plhs[4]);
  
  
  /* setup interface interpreter and wrap around "x" another structure */
   MATLABSetupInterpreter( &ii );
   xx = mv_MultiVectorWrap( &ii, x, 0 /* "0" means doesn't own data */);
  
/* set pointers to lapack functions */
   blap_fn.dpotrf = dpotrf_;
   blap_fn.dsygv = dsygv_;
   
   /* handle constraints */
   if (!mxIsEmpty(prhs[4]))
   {
     y = MatlabMultiVectorCreate(prhs[4]);
     yy = mv_MultiVectorWrap( &ii, y, 0 /* "0" means doesn't own data */);
   }
   
/* execute lobpcg */
   exitflag = lobpcg_solve_double( xx,
	             (void*)prhs[1],      /* operatorA */
                 MatMultiVec,  
	             (void*)prhs[2],      /* operatorB */
	             mxIsEmpty(prhs[2])? NULL : MatMultiVec,
	             (void*)prhs[3],      /* operatorT */
	             mxIsEmpty(prhs[3])? NULL : MatMultiVec,
	             mxIsEmpty(prhs[4])? NULL : yy,   /* constraints */
                 blap_fn,
	             lobpcg_tol,
	             maxit,
	             verbosity_level,
	             &iterations,

/* eigenvalues; "lambda_values" should point to array  containing <blocksize> doubles where <blocksi
ze> is the width of multivector "blockVectorX" */
                 eigs,

/* eigenvalues history; a pointer to the entries of the  <blocksize>-by-(<maxIterations>+1) matrix s
tored
in  fortran-style. (i.e. column-wise) The matrix may be  a submatrix of a larger matrix, see next
argument; If you don't need eigenvalues history, provide NULL in this entry */
                 eigs_hist,

/* global height of the matrix (stored in fotran-style)  specified by previous argument */
                 width,

/* residual norms; argument should point to array of <blocksize> doubles */
                 resid,

/* residual norms history; a pointer to the entries of the  <blocksize>-by-(<maxIterations>+1) matri
x
stored in  fortran-style. (i.e. column-wise) The matrix may be  a submatrix of a larger matrix, see
next
argument If you don't need residual norms history, provide NULL in this entry */
                 resid_hist,

/* global height of the matrix (stored in fotran-style)  specified by previous argument */
                 width
	);
   
   plhs[2] = mxCreateLogicalScalar(exitflag!=0);
   plhs[5] = mxCreateDoubleScalar(iterations);
   /* extract blockVector with calculated eigenvectors */
   plhs[0] = MatlabMultiVectorGetCopyOfBlockVector(x);
   
   /* clean up */
   MatlabMultiVectorDestroy(x);
   mv_MultiVectorDestroy(xx);
   mxFree(resid);
   if (!mxIsEmpty(prhs[4]))
   {
     MatlabMultiVectorDestroy(y);
     mv_MultiVectorDestroy(yy);
   }
}
