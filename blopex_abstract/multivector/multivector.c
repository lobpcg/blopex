/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ University of Colorado at Denver                                  */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "multivector.h"

/* -----------------------------------------------------
   abstract multivector
   ----------------------------------------------------- */
struct mv_MultiVector
{
  void* data;      /* the pointer to the actual multivector */
  int   ownsData;

  mv_InterfaceInterpreter* interpreter; /* a structure that defines
                                           multivector operations */
} ;
/* ------------------------------------------------------
   mv_MultiVectorGetData                          generic
   ------------------------------------------------------ */
void *
mv_MultiVectorGetData (mv_MultiVectorPtr x)
{
  assert (x!=NULL);
  return x->data;
}
/* -------------------------------------------------------
   mv_MultiVectorWrap                             generic
   ------------------------------------------------------- */
mv_MultiVectorPtr
mv_MultiVectorWrap( mv_InterfaceInterpreter* ii, void * data, int ownsData )
{
  mv_MultiVectorPtr x;

  x = (mv_MultiVectorPtr) malloc(sizeof(struct mv_MultiVector));
  assert( x != NULL );

  x->interpreter = ii;
  x->data = data;
  x->ownsData = ownsData;

  return x;
}
/* ---------------------------------------------------------------
   mv_MultiVectorCreateFromSampleVector      not used      generic
   --------------------------------------------------------------- */
mv_MultiVectorPtr
mv_MultiVectorCreateFromSampleVector( void* ii_, int n, void* sample ) {

  mv_MultiVectorPtr x;
  mv_InterfaceInterpreter* ii = (mv_InterfaceInterpreter*)ii_;

  x = (mv_MultiVectorPtr) malloc(sizeof(struct mv_MultiVector));
  assert( x != NULL );

  x->interpreter = ii;
  x->data = (ii->CreateMultiVector)( ii, n, sample );
  x->ownsData = 1;

  return x;
}
/* ------------------------------------------------------------------
   mv_MultiVectorCreateCopy                                   generic
   ------------------------------------------------------------------ */
mv_MultiVectorPtr
mv_MultiVectorCreateCopy( mv_MultiVectorPtr x, int copyValues ) {

  mv_MultiVectorPtr y;
  void* data;
  mv_InterfaceInterpreter* ii;

  assert( x != NULL );
  ii = x->interpreter;

  y = (mv_MultiVectorPtr) malloc(sizeof(struct mv_MultiVector));
  assert( y != NULL );

  data = (ii->CopyCreateMultiVector)( x->data, copyValues );

  y->interpreter = ii;
  y->data = data;
  y->ownsData = 1;

  return y;
}
/* ------------------------------------------------------------------
   mv_MultiVectorDestroy                                      generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorDestroy( mv_MultiVectorPtr v) {

  if ( v == NULL )
    return;

  if ( v->ownsData )
    (v->interpreter->DestroyMultiVector)( v->data );
  free( v );
}
/* ------------------------------------------------------------------
   mv_MultiVectorSetMask                                      generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorSetMask( mv_MultiVectorPtr v, int* mask ) {

  assert( v != NULL );
  (v->interpreter->SetMask)( v->data, mask );
}
/* ------------------------------------------------------------------
   mv_MultiVectorWidth                                        generic
   ------------------------------------------------------------------ */
int
mv_MultiVectorWidth( mv_MultiVectorPtr v ) {

  if ( v == NULL )
    return 0;

  return (v->interpreter->Width)( v->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorHeight                 not used              generic
   ------------------------------------------------------------------ */
int
mv_MultiVectorHeight( mv_MultiVectorPtr v ) {

  if ( v == NULL )
    return 0;

  return (v->interpreter->Height)(v->data);
}
/* ------------------------------------------------------------------
   mv_MultiVectorClear                                        generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorClear( mv_MultiVectorPtr v ) {

  assert( v != NULL );
  (v->interpreter->ClearMultiVector)( v->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorSetRandom                                    generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorSetRandom( mv_MultiVectorPtr v, int seed ) {

  assert( v != NULL );
  (v->interpreter->SetRandomVectors)( v->data, seed );
}
/* ------------------------------------------------------------------
   mv_MultiVectorCopy                                         generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorCopy( mv_MultiVectorPtr src, mv_MultiVectorPtr dest ) {

  assert( src != NULL && dest != NULL );
  (src->interpreter->CopyMultiVector)( src->data, dest->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorAxpy        y=y+a*x                          generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorAxpy( double a, mv_MultiVectorPtr x, mv_MultiVectorPtr y ) {

  assert( x != NULL && y != NULL );
  (x->interpreter->MultiAxpy)( a, x->data, y->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorByMultiVector    xy=x'*y                     generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorByMultiVector( mv_MultiVectorPtr x,
                             mv_MultiVectorPtr y,
                             int xyGHeight, int xyHeight,
                             int xyWidth, void* xy ) {

  assert( x != NULL && y != NULL );
  (x->interpreter->MultiInnerProd)
    ( x->data, y->data, xyGHeight, xyHeight, xyWidth, xy );
}
/* ------------------------------------------------------------------
   mv_MultiVectorByMultiVectorDiag     d=diag(x'*y)           generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorByMultiVectorDiag( mv_MultiVectorPtr x,
                                 mv_MultiVectorPtr y,
                                 int* mask, int n, void* d ) {

  assert( x != NULL && y != NULL );
  (x->interpreter->MultiInnerProdDiag)( x->data, y->data, mask, n, d );
}
/* ------------------------------------------------------------------
   mv_MultiVectorByMatrix              y=x*r                  generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorByMatrix( mv_MultiVectorPtr x,
                        int rGHeight, int rHeight, int rWidth,
                        void* rVal,
                        mv_MultiVectorPtr y ) {

  assert( x != NULL && y != NULL );
  (x->interpreter->MultiVecMat)
    ( x->data, rGHeight, rHeight, rWidth, rVal, y->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorXapy                y=y+x*a      not used    generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorXapy( mv_MultiVectorPtr x,
                    int rGHeight, int rHeight, int rWidth,
                    void* rVal,
                    mv_MultiVectorPtr y ) {

  assert( x != NULL && y != NULL );
  (x->interpreter->MultiXapy)
    ( x->data, rGHeight, rHeight, rWidth, rVal, y->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorByDiagonal            y=x*d(mask)            generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorByDiagonal( mv_MultiVectorPtr x,
                          int* mask, int n, void* d,
                          mv_MultiVectorPtr y ) {

  /* y = x*d */

  assert( x != NULL && y != NULL );
  (x->interpreter->MultiVecMatDiag)( x->data, mask, n, d, y->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorEval                     not used            generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorEval( void (*f)( void*, void*, void* ), void* par,
               mv_MultiVectorPtr x, mv_MultiVectorPtr y ) {

  /* y = f(x) computed vector-wise */

  assert( x != NULL && y != NULL );
  (x->interpreter->Eval)( f, par, x->data, y->data );
}
/* ------------------------------------------------------------------
   mv_MultiVectorPrint                                        generic
   ------------------------------------------------------------------ */
void
mv_MultiVectorPrint( mv_MultiVectorPtr x,char * tag,int limit ) {

  assert( x != NULL );
  (x->interpreter->MultiPrint)( x->data, tag, limit );
}
