/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.2) MIT / Apache-2.0 License                   */
/* @@@ Copyright 2010-2019 Team https://github.com/lobpcg/blopex       */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#include "mex.h"
#include "matrix.h"

#include "matlab_interface.h"
#include "multivector_for_matlab.h"

#include "fortran_matrix.h"
#include "lobpcg.h"

#define REQUIRED_NUM_INPUTS    8
#define REQUIRED_NUM_OUTPUTS   6

int BLOPEX_MATLAB_COMPLEX = 0;

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
/* define complex variables */
  komplex*                     eigs_complex;
  komplex*                     eigs_hist_complex;
  int                          k;
/* define double variables */
  double*                      eigs;
  double*                      eigs_hist;

  double*                      resid;
  double*                      resid_hist;
  int                          maxit;
  int                          verbosity_level;
  BlopexInt                    width;
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
 
/* determine if for complex arithematic */
  if (mxIsComplex(prhs[1]))
    BLOPEX_MATLAB_COMPLEX = 1;
  else
    BLOPEX_MATLAB_COMPLEX = 0;

  printf("Complex=%d\n",BLOPEX_MATLAB_COMPLEX);
 
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
if (BLOPEX_MATLAB_COMPLEX) 
  eigs_complex = (komplex*)malloc(sizeof(komplex)*width);
else {
  plhs[1] = mxCreateDoubleMatrix(width, 1, mxREAL);
  eigs = mxGetPr(plhs[1]); 
 } 
   /* request mem. for eig-history */
if (BLOPEX_MATLAB_COMPLEX)
  eigs_hist_complex = (komplex*)malloc(sizeof(komplex)*width*(maxit+1));
else {
  plhs[3] = mxCreateDoubleMatrix(width, maxit+1, mxREAL);
  eigs_hist = mxGetPr(plhs[3]);
}
  
   /* request memory for resid. norms */
  resid = mxCalloc(width, sizeof(double));
  
  /* req. mem. for resid-hist */
  plhs[4] = mxCreateDoubleMatrix(width, maxit+1, mxREAL);
  resid_hist = mxGetPr(plhs[4]);
  
  
  /* setup interface interpreter and wrap around "x" another structure */
   MATLABSetupInterpreter( &ii );
   xx = mv_MultiVectorWrap( &ii, x, 0 /* "0" means doesn't own data */);
  
/* set pointers to lapack functions */
if (BLOPEX_MATLAB_COMPLEX) {
   blap_fn.zpotrf = zpotrf_redirect;
   blap_fn.zhegv = zhegv_redirect;
}
else {
   blap_fn.dpotrf = dpotrf_redirect;
   blap_fn.dsygv = dsygv_redirect;
}
   
   /* handle constraints */
   if (!mxIsEmpty(prhs[4]))
   {
     y = MatlabMultiVectorCreate(prhs[4]);
     yy = mv_MultiVectorWrap( &ii, y, 0 /* "0" means doesn't own data */);
   }
   
/* execute lobpcg */
if (BLOPEX_MATLAB_COMPLEX) {
   printf("Executing complex");
   exitflag = lobpcg_solve_complex( xx,
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
                 eigs_complex,

/* eigenvalues history; a pointer to the entries of the  <blocksize>-by-(<maxIterations>+1) matrix s
tored
in  fortran-style. (i.e. column-wise) The matrix may be  a submatrix of a larger matrix, see next
argument; If you don't need eigenvalues history, provide NULL in this entry */
                 eigs_hist_complex,

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
}
else {
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
}
   
   plhs[2] = mxCreateLogicalScalar(exitflag!=0);
   plhs[5] = mxCreateDoubleScalar(iterations);
   /* extract blockVector with calculated eigenvectors */
   plhs[0] = MatlabMultiVectorGetCopyOfBlockVector(x);
   /* extract computed eigenvalues(real) for complex arithmetic case */
if (BLOPEX_MATLAB_COMPLEX) {   
   plhs[1] = mxCreateDoubleMatrix(width, 1, mxREAL);
   for (k=0; k<width; k++)
     mxGetPr(plhs[1])[k] = eigs_complex[k].real;
   free(eigs_complex);
} 
   /* extract eigenvalue history for complex arithmetic case */
if (BLOPEX_MATLAB_COMPLEX) {
  plhs[3] = mxCreateDoubleMatrix(width, maxit+1, mxCOMPLEX);
  for (k=0; k<width*(maxit+1); k++)
  {
    mxGetPr(plhs[3])[k] = eigs_hist_complex[k].real;
    mxGetPi(plhs[3])[k] = eigs_hist_complex[k].imag;
  }  
  free(eigs_hist_complex); 
}    
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
