/*BHEADER**********************************************************************
 * Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * This file is part of HYPRE.  See file COPYRIGHT for details.
 *
 * HYPRE is free software; you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * $Revision: 2.8 $
 ***********************************************************************EHEADER*/




#include "_hypre_struct_ls.h"
#include "temp_multivector.h"

int 
hypre_StructVectorSetRandomValues( hypre_StructVector *vector,
                                   int seed )
{
   int    ierr = 0;

   hypre_Box          *v_data_box;
                    
   int                 vi;
   double             *vp;

   hypre_BoxArray     *boxes;
   hypre_Box          *box;
   hypre_Index         loop_size;
   hypre_IndexRef      start;
   hypre_Index         unit_stride;

   int                 i;
   int                 loopi, loopj, loopk;

   /*-----------------------------------------------------------------------
    * Set the vector coefficients
    *-----------------------------------------------------------------------*/

   srand( seed );

   hypre_SetIndex(unit_stride, 1, 1, 1);
 
   boxes = hypre_StructGridBoxes(hypre_StructVectorGrid(vector));
   hypre_ForBoxI(i, boxes)
      {
         box      = hypre_BoxArrayBox(boxes, i);
         start = hypre_BoxIMin(box);

         v_data_box =
            hypre_BoxArrayBox(hypre_StructVectorDataSpace(vector), i);
         vp = hypre_StructVectorBoxData(vector, i);
 
         hypre_BoxGetSize(box, loop_size);

         hypre_BoxLoop1Begin(loop_size,
                             v_data_box, start, unit_stride, vi);
#define HYPRE_BOX_SMP_PRIVATE loopk,loopi,loopj,vi 
#include "hypre_box_smp_forloop.h"
         hypre_BoxLoop1For(loopi, loopj, loopk, vi)
            {
               vp[vi] = 2.0*rand()/RAND_MAX - 1.0;
            }
         hypre_BoxLoop1End(vi);
      }

   return ierr;
}

int
hypre_StructSetRandomValues( void* v, int seed ) {

  return hypre_StructVectorSetRandomValues( (hypre_StructVector*)v, seed );
}

/* these WRAPPER functions added for blopex complex additions, 
 * for pointers that were double but are now void to 
 * support both double and complex
 */

int hypre_StructKrylovInnerProd_WRAPPER  ( void *x, void *y, void *result )
{
   *((double *)result) = hypre_StructKrylovInnerProd (x,y);
   return 0;
}

int hypre_StructKrylovAxpy_WRAPPER ( void * alpha, void *x, void *y )
{
   hypre_StructKrylovAxpy( *((double *)alpha), x, y );
   return 0;
}

int
HYPRE_StructSetupInterpreter( mv_InterfaceInterpreter *i )
{
  i->CreateVector = hypre_StructKrylovCreateVector;
  i->DestroyVector = hypre_StructKrylovDestroyVector; 
  i->InnerProd = hypre_StructKrylovInnerProd_WRAPPER; 
  i->CopyVector = hypre_StructKrylovCopyVector;
  i->ClearVector = hypre_StructKrylovClearVector;
  i->SetRandomValues = hypre_StructSetRandomValues;
  i->ScaleVector = hypre_StructKrylovScaleVector;
  i->Axpy = hypre_StructKrylovAxpy_WRAPPER;

  i->CreateMultiVector = mv_TempMultiVectorCreateFromSampleVector;
  i->CopyCreateMultiVector = mv_TempMultiVectorCreateCopy;
  i->DestroyMultiVector = mv_TempMultiVectorDestroy;

  i->Width = mv_TempMultiVectorWidth;
  i->Height = mv_TempMultiVectorHeight;
  i->SetMask = mv_TempMultiVectorSetMask;
  i->CopyMultiVector = mv_TempMultiVectorCopy;
  i->ClearMultiVector = mv_TempMultiVectorClear;
  i->SetRandomVectors = mv_TempMultiVectorSetRandom;
  i->MultiInnerProd = mv_TempMultiVectorByMultiVector;
  i->MultiInnerProdDiag = mv_TempMultiVectorByMultiVectorDiag;
  i->MultiVecMat = mv_TempMultiVectorByMatrix;
  i->MultiVecMatDiag = mv_TempMultiVectorByDiagonal;
  i->MultiAxpy = mv_TempMultiVectorAxpy;
  i->MultiXapy = mv_TempMultiVectorXapy;
  i->Eval = mv_TempMultiVectorEval;

  return 0;
}

int
HYPRE_StructSetupMatvec(HYPRE_MatvecFunctions * mv)
{
  mv->MatvecCreate = hypre_StructKrylovMatvecCreate;
  mv->Matvec = hypre_StructKrylovMatvec;
  mv->MatvecDestroy = hypre_StructKrylovMatvecDestroy;

  mv->MatMultiVecCreate = NULL;
  mv->MatMultiVec = NULL;
  mv->MatMultiVecDestroy = NULL;

  return 0;
}
