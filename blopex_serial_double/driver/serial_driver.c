/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */
/* test driver for double serial implementation of lobpcg */

#define MAXIT 100
#define TOL 1e-6

#include "fortran_matrix.h"
#include "fortran_interpreter.h"
#include "lobpcg.h"
#include "multi_vector.h"
#include "multivector.h"
#include "pcg_multi.h"
#include "interpreter.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int dsygv_ (int *itype, char *jobz, char *uplo, int *
                    n, double *a, int *lda, double *b, int *ldb,
                    double *w, double *work, int *lwork, int *info);

int dpotrf_ (char *uplo, int *n, double *a, int *
                    lda, int *info);

int main(void)
{
   serial_Multi_Vector * x;
   serial_Multi_Vector * operatorA;
   mv_MultiVectorPtr xx;
   double * eigs;
   double * resid;
   int iterations;
   lobpcg_Tolerance lobpcg_tol;
   mv_InterfaceInterpreter ii;
   lobpcg_BLASLAPACKFunctions blap_fn;

   int MV_HEIGHT, MV_WIDTH;

/* create operatorA */
   operatorA = serial_Multi_VectorCreate(40,40);
   serial_Multi_VectorInitialize(operatorA);
   double kzero = 0.0;
   serial_Multi_VectorSetConstantValues( operatorA, kzero);
   double * pc;
   int di;
   pc = (double *)operatorA->data;
   for (di=0;di<40;di++) {
     *pc = di+1;
     pc=pc+41;
   }
   serial_Multi_VectorPrint(operatorA,"operatorA",4);

   MV_HEIGHT = operatorA->size;
   MV_WIDTH =  5;

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

/* execute lobpcg to solve for Ax=(lambda)Bx
   with initial guess of e-vectors xx
   and preconditioner T
   number of vectors in xx determines number of e-values eig to solve for
   solving for the smallest e-values under constraints Y
   execution stops when e-values with tolerances or after max iterations */

   printf("Call lobpcg solve double\n");

    lobpcg_solve_double(
              xx,           /*input-initial guess of e-vectors */
   (void *) operatorA,      /*input-matrix A */
          MatMultiVec,      /*input-operator A */
          NULL,             /*input-matrix B */
          NULL,             /*input-operator B */
          NULL,             /*input-matrix T */
          NULL,             /*input-operator T */
          NULL,             /*input-matrix Y */
              blap_fn,      /*input-lapack functions */
              lobpcg_tol,   /*input-tolerances */
              MAXIT,        /*input-max iterations */
              2,            /*input-verbosity level */

              &iterations,  /*output-actual iterations */
              eigs,         /*output-eigenvalues */
          NULL,             /*output-eigenvalues history */
              0,            /*output-history global height */
              resid,        /*output-residual norms */
           NULL,            /*output-residual norms history */
              0             /*output-history global height  */
    );

/* print eigenvectors */
/* serial_Multi_VectorPrint(x, "eigenvectors", 0);  */

/* destroy multivector and other objects */

   serial_Multi_VectorDestroy(x);
   mv_MultiVectorDestroy(xx);
   free(eigs);
   free(resid);

   return 0;
}
