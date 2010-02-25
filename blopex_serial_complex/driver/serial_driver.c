/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 2.0) LGPL Version 3 or above.  See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team http://code.google.com/p/blopex/     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov, Don McCuan */
/* test driver for complex implementation of lobpcg */

#define MAXIT 1000

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
#include <time.h>

   BlopexInt zpotrf_ (char *uplo, BlopexInt *n,
                komplex *a, BlopexInt *lda, BlopexInt *info);
   BlopexInt zhegv_  (BlopexInt *itype, char *jobz, char *uplo, BlopexInt * n,
                komplex *a, BlopexInt *lda, komplex *b, BlopexInt *ldb,
                double *w, komplex *work, BlopexInt *lwork,
                double * rwork,BlopexInt *info);

BlopexInt main(BlopexInt argc,char *argv[])
{

   serial_Multi_Vector * x;
   serial_Multi_Vector * operatorA;
   mv_MultiVectorPtr xx;
   komplex * eigs;
   double  * resid;
   int iterations;
   lobpcg_Tolerance lobpcg_tol;
   mv_InterfaceInterpreter ii;
   lobpcg_BLASLAPACKFunctions blap_fn;

   BlopexInt MV_HEIGHT, MV_WIDTH;

   if (argc != 4) {
     printf("Must have 3 arguments. Try:\n");
     printf("./serial_driver L test1a.txt test1x.txt\n");
     printf(" or \n");
     printf("./serial_driver R 100 10\n");
     return(1);
   }

   /* ====== argv = L Aname Xname   */
   if (*argv[1]=='L') {
     printf("Vector Load into operatorA\n");
     operatorA = serial_Multi_VectorLoad(argv[2]);
     serial_Multi_VectorPrint(operatorA,"operatorA",5);

     MV_HEIGHT = operatorA->size;
     printf("Vector Load into x");
     x = serial_Multi_VectorLoad(argv[3]);
     serial_Multi_VectorPrint(x,"x",5);
     MV_WIDTH = x->num_vectors;
   }
   /* ===== argv = R height width   */
   else {
     printf("Create random vectors\n");
     MV_HEIGHT = atoi(argv[2]);
     MV_WIDTH  = atoi(argv[3]);

     operatorA = serial_Multi_VectorCreate(MV_HEIGHT,MV_HEIGHT);
     serial_Multi_VectorInitialize(operatorA);
     serial_Multi_VectorSetRandomValues(operatorA, 1);
     serial_Multi_VectorPrint(operatorA,"operatorA",5);
     serial_Multi_VectorSymmetrize(operatorA);
     serial_Multi_VectorPrint(operatorA,"operatorA",5);
/*
     komplex kzero = {0.0, 0.0};
     serial_Multi_VectorSetConstantValues( operatorA, kzero);
     komplex * pc;
     BlopexInt di;
     pc = (komplex *)operatorA->data;
     for (di=0;di<MV_HEIGHT;di++) {
       pc->real = di+1;
       pc=pc+MV_HEIGHT+1;
*/
     x = serial_Multi_VectorCreate(MV_HEIGHT, MV_WIDTH);
     serial_Multi_VectorInitialize(x);
     serial_Multi_VectorSetRandomValues(x, 1);
     serial_Multi_VectorPrint(x,"x",5);
   }

/* get memory for eigenvalues, eigenvalue history, residual norms, residual norms history */
   eigs = (komplex *)malloc(sizeof(komplex)*MV_WIDTH);
   resid = (double *)malloc(sizeof(double)*MV_WIDTH);

/* set tolerances */
   lobpcg_tol.absolute = 1e-6;
   lobpcg_tol.relative = 1e-50;

/* setup interface interpreter and wrap around "x" another structure */
   SerialSetupInterpreter( &ii );
   xx = mv_MultiVectorWrap( &ii, x, 0);

/* set pointers to lapack functions */
   blap_fn.zpotrf = zpotrf_;
   blap_fn.zhegv = zhegv_;

/* execute lobpcg to solve for Ax=(lambda)Bx
   with initial guess of e-vectors xx
   and preconditioner T
   number of vectors in xx determines number of e-values eig to solve for
   solving for the smallest e-values under constraints Y
   execution stops when e-values with tolerances or after max iterations */

   printf("Call lobpcg solve complex\n");

   clock_t start,end;
   double cpu_time_used;

   start = clock();

   lobpcg_solve_complex(
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
              1,            /*input-verbosity level */

              &iterations,  /*output-actual iterations */
              eigs,         /*output-eigenvalues */
          NULL,             /*output-eigenvalues history */
              0,            /*output-history global height */
              resid,        /*output-residual norms */
           NULL,            /*output-residual norms history */
              0             /*output-history global height  */
    );

   end = clock();
   cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

   printf("CPU time used = %22.16e \n",cpu_time_used);

/* print eigenvector */
/*
   serial_Multi_VectorPrint(x,"eigenvectors",0);
*/

/* destroy multivector and other objects */
   serial_Multi_VectorDestroy(x);
   mv_MultiVectorDestroy(xx);
   free(eigs);
   free(resid);

   return 0;
}
