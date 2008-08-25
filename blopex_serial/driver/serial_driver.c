/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */
/* test driver for serial implementation of lobpcg */

#define MV_HEIGHT 100
#define MV_WIDTH 10
#define MAXIT 100
#define TOL 1e-6

#include "lobpcg.h"
#include "multi_vector.h"
#include "multivector.h"
#include "pcg_multi.h"
#include "interpreter.h"

#include <stdlib.h>

int dsygv_ (int *itype, char *jobz, char *uplo, int *
                    n, double *a, int *lda, double *b, int *ldb,
                    double *w, double *work, int *lwork, int *info);
                    
int dpotrf_ (char *uplo, int *n, double *a, int *
                    lda, int *info);

int main(void)
{
   serial_Multi_Vector * x;
   mv_MultiVectorPtr xx; 
   double * eigs;
   double * resid;
   int iterations;
   lobpcg_Tolerance lobpcg_tol;
   mv_InterfaceInterpreter ii;
   lobpcg_BLASLAPACKFunctions blap_fn;

/* create multivector */
   x = serial_Multi_VectorCreate(MV_HEIGHT, MV_WIDTH);
   serial_Multi_VectorInitialize(x);
      
/* fill it with random numbers */
   serial_Multi_VectorSetRandomValues(x, 1);

/* get memory for eigenvalues, eigenvalue history, residual norms, residual norms history */

  /* request memory for eig-vals */
   eigs = (double *)malloc(sizeof(double)*MV_WIDTH);
   
   /* request memory for resid. norms */
   resid = (double *)malloc(sizeof(double)*MV_WIDTH);
   
/* set tolerances */
   lobpcg_tol.absolute = TOL;
   lobpcg_tol.relative = 1e-50;
   
/* setup interface interpreter and wrap around "x" another structure */

   SerialSetupInterpreter( &ii );
   xx = mv_MultiVectorWrap( &ii, x, 0);
  
/* set pointers to lapack functions */
   blap_fn.dpotrf = dpotrf_;
   blap_fn.dsygv = dsygv_;
  
/* execute lobpcg */
   lobpcg_solve( xx,
	         NULL,
              MatMultiVec,
	      NULL,
	      NULL,
	      NULL,
	      NULL,
	      NULL,
              blap_fn,
	      lobpcg_tol,
	      MAXIT,
	      1,
	      &iterations,

/* eigenvalues; "lambda_values" should point to array  containing <blocksize> doubles where <blocksi
ze> is the width of multivector "blockVectorX" */
              eigs,

/* eigenvalues history; a pointer to the entries of the  <blocksize>-by-(<maxIterations>+1) matrix s
tored
in  fortran-style. (i.e. column-wise) The matrix may be  a submatrix of a larger matrix, see next
argument; If you don't need eigenvalues history, provide NULL in this entry */
              NULL,

/* global height of the matrix (stored in fotran-style)  specified by previous argument */
              0,

/* residual norms; argument should point to array of <blocksize> doubles */
              resid,

/* residual norms history; a pointer to the entries of the  <blocksize>-by-(<maxIterations>+1) matri
x
stored in  fortran-style. (i.e. column-wise) The matrix may be  a submatrix of a larger matrix, see
next
argument If you don't need residual norms history, provide NULL in this entry */
              NULL,

/* global height of the matrix (stored in fotran-style)  specified by previous argument */
              0
	);

/* print eigenvectors to file */
   serial_Multi_VectorPrint(x, "eigenvectors");
   
/* destroy multivector and other objects */
   
   serial_Multi_VectorDestroy(x);
   mv_MultiVectorDestroy(xx);
   free(eigs); 
   free(resid); 
   
   return 0;
}
