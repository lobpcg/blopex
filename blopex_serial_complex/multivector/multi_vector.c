/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1)                                              */
/* @@@ University of Colorado at Denver                                  */
/* @@@ Copyright 2008 Merico Argentati, Andrew Knyazev,                  */
/* @@@ Ilya Lashuk, Evgueni Ovtchinnikov, and Don McCuan                 */
/* @@@ LGPL Version 3 or above.  See www.gnu.org.                        */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov, Don McCuan */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

typedef struct {double real, imag;} komplex;

void complex_multiply(komplex* A, komplex* B, komplex* C);
void complex_add(komplex* A, komplex* B, komplex* C);
void complex_subtract(komplex* A, komplex* B, komplex* C);
void complex_divide(komplex* A, komplex* B, komplex* C);

#include "multi_vector.h"

/*--------------------------------------------------------------------------
 * serial_Multi_VectorCreate                                         generic
 *--------------------------------------------------------------------------*/

serial_Multi_Vector *
serial_Multi_VectorCreate( int size, int num_vectors  )
{
   serial_Multi_Vector *mvector;

   mvector = (serial_Multi_Vector *) malloc (sizeof(serial_Multi_Vector));

   serial_Multi_VectorNumVectors(mvector) = num_vectors;
   serial_Multi_VectorSize(mvector) = size;

   serial_Multi_VectorOwnsData(mvector) = 1;
   serial_Multi_VectorData(mvector) = NULL;


   mvector->num_active_vectors=0;
   mvector->active_indices=NULL;

   return mvector;
}
/*--------------------------------------------------------------------------
 * serial_Multi_VectorInitialize                                    complex
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorInitialize( serial_Multi_Vector *mvector )
{
   int    ierr = 0, i, size, num_vectors;

   size        = serial_Multi_VectorSize(mvector);
   num_vectors = serial_Multi_VectorNumVectors(mvector);

   if (NULL==serial_Multi_VectorData(mvector))
      serial_Multi_VectorData(mvector) = malloc (sizeof(komplex)*size*num_vectors);

   /* now we create a "mask" of "active" vectors; initially all vectors are active */
   if (NULL==mvector->active_indices)
    {
         mvector->active_indices=(int *) malloc(sizeof(int)*num_vectors);

         for (i=0; i<num_vectors; i++)
            mvector->active_indices[i]=i;

         mvector->num_active_vectors=num_vectors;
    }

   return ierr;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorSetDataOwner                                   generic
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorSetDataOwner( serial_Multi_Vector *mvector, int owns_data )
{
   int    ierr=0;

   serial_Multi_VectorOwnsData(mvector) = owns_data;

   return ierr;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorDestroy                                        generic
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorDestroy( serial_Multi_Vector *mvector )
{
   int    ierr=0;

   if (NULL!=mvector)
   {
      if (serial_Multi_VectorOwnsData(mvector) && NULL!=serial_Multi_VectorData(mvector))
         free( serial_Multi_VectorData(mvector) );

      if (NULL!=mvector->active_indices)
            free(mvector->active_indices);

      free(mvector);
   }
   return ierr;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorSetMask                                        generic
 *--------------------------------------------------------------------------*/
 int
 serial_Multi_VectorSetMask(serial_Multi_Vector *mvector, int * mask)
 {
   /* this routine accepts mask in "zeros and ones format, and converts it to the one used in
   the structure "serial_Multi_Vector" */
   int  num_vectors = mvector->num_vectors;
   int i;


   /* may be it's better to just check if it is not null, and throw an error, if it is? */
   if (mvector->active_indices==NULL)
      mvector->active_indices=(int *) malloc(sizeof(int)*num_vectors);

   mvector->num_active_vectors=0;

   if (mask!=NULL)
      for (i=0; i<num_vectors; i++)
      {
         if ( mask[i] )
            mvector->active_indices[mvector->num_active_vectors++]=i;
      }
   else
      for (i=0; i<num_vectors; i++)
         mvector->active_indices[mvector->num_active_vectors++]=i;

   return 0;
 }


/*--------------------------------------------------------------------------
 * serial_Multi_VectorSetConstantValues                             complex
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorSetConstantValues( serial_Multi_Vector *v,
                                       komplex value)
{
   komplex  *vector_data = (komplex *) serial_Multi_VectorData(v);
   int      size        = serial_Multi_VectorSize(v);
   int      i, j, start_offset, end_offset;

   for (i = 0; i < v->num_active_vectors; i++)
   {
      start_offset = v->active_indices[i]*size;
      end_offset = start_offset+size;

      for (j=start_offset; j < end_offset; j++)
         vector_data[j]= value;
   }

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorSetRandomValues                               complex
 *
 *     returns vector of values randomly distributed between -1.0 and +1.0
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorSetRandomValues( serial_Multi_Vector *v, int seed)
{
   komplex  *vector_data = (komplex *) serial_Multi_VectorData(v);
   int      size        = serial_Multi_VectorSize(v);
   int      i, j, start_offset, end_offset;

   srand48(seed);
   printf("random size %d\n",size);
   printf("random active vectors %d\n",v->num_active_vectors);
   for (i = 0; i < v->num_active_vectors; i++)
   {
      start_offset = v->active_indices[i]*size;
      end_offset = start_offset+size;

      for (j=start_offset; j < end_offset; j++) {
         vector_data[j].real= 2.0 * drand48() - 1.0;
         vector_data[j].imag= 2.0 * drand48() - 1.0;;
      }
   }

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorCopy   y=x  using indices                     complex
 *
 * y should have already been initialized at the same size as x
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorCopy( serial_Multi_Vector *x, serial_Multi_Vector *y)
{
   komplex *x_data;
   komplex *y_data;
   int i;
   int size;
   int num_bytes;
   int num_active_vectors;
   komplex * dest;
   komplex * src;
   int * x_active_ind;
   int * y_active_ind;

   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);

   num_active_vectors = x->num_active_vectors;
   size = x->size;
   num_bytes = size*sizeof(komplex);
   x_data = (komplex *) x->data;
   y_data = (komplex *) y->data;
   x_active_ind=x->active_indices;
   y_active_ind=y->active_indices;

   for (i=0; i<num_active_vectors; i++)
   {
      src=x_data + size * x_active_ind[i];
      dest = y_data + size * y_active_ind[i];

      memcpy(dest,src,num_bytes);
   }

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorCopyWithoutMask              y=x              complex
 * copies data from x to y without using indices
 * y should have already been initialized at the same size as x
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorCopyWithoutMask(serial_Multi_Vector *x , serial_Multi_Vector *y)
{
   int byte_count;

   assert (x->size == y->size && x->num_vectors == y->num_vectors);

   byte_count = sizeof(komplex) * x->size * x->num_vectors;

/* threading not done here since it's not done (reason?) in serial_VectorCopy
      from vector.c */

   memcpy(y->data,x->data,byte_count);

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorAxpy      y=alpha*x+y using indices           complex
 *
 * call seq < MultiVectorAxpy < i->MultiAxpy < mv_MultiVectorAxpy
 * alpha is always 1 or -1 double
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorAxpy( double           alpha,
                          serial_Multi_Vector *x,
                          serial_Multi_Vector *y)
{
   komplex  *x_data;
   komplex  *y_data;
   komplex * src;
   komplex * dest;
   int * x_active_ind;
   int * y_active_ind;
   int i, j;
   int size;
   int num_active_vectors;

   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);

   x_data = (komplex *) x->data;
   y_data = (komplex *) y->data;
   size = x->size;
   num_active_vectors = x->num_active_vectors;
   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;

   for(i=0; i<num_active_vectors; i++)
   {
      src = x_data + x_active_ind[i]*size;
      dest = y_data + y_active_ind[i]*size;

      for (j=0; j<size; j++) {
         dest[j].real += alpha*src[j].real;
         dest[j].imag += alpha*src[j].imag;
      }
   }

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorByDiag:    y=x*alpha using indices            complex
 *
 * if y and x are mxn then alpha is nx1
 * call seq < MultiVectorByDiagonal < i->MultiVecMatDiag < mv_MultiVectorByDiagonal
 *--------------------------------------------------------------------------*/
int
serial_Multi_VectorByDiag( serial_Multi_Vector *x,
                            int                 *mask,
                            int                 n,
                            komplex             *alpha,
                            serial_Multi_Vector *y)
{
   komplex  *x_data;
   komplex  *y_data;
   int      size;
   int      num_active_vectors;
   int      i,j;
   komplex  *dest;
   komplex  *src;
   int * x_active_ind;
   int * y_active_ind;
   int * al_active_ind;
   int num_active_als;
   komplex * current_alpha;

   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);

   /* build list of active indices in alpha */

   al_active_ind = (int *) malloc(sizeof(int)*n);
   num_active_als = 0;

   if (mask!=NULL)
      for (i=0; i<n; i++)
      {
         if (mask[i])
            al_active_ind[num_active_als++]=i;
      }
   else
      for (i=0; i<n; i++)
         al_active_ind[num_active_als++]=i;

   assert (num_active_als==x->num_active_vectors);

   x_data = (komplex *) x->data;
   y_data = (komplex *) y->data;
   size = x->size;
   num_active_vectors = x->num_active_vectors;
   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;

   for(i=0; i<num_active_vectors; i++)
   {
      src = x_data + x_active_ind[i]*size;
      dest = y_data + y_active_ind[i]*size;
      current_alpha=&alpha[ al_active_ind[i] ];

      for (j=0; j<size; j++)
         complex_multiply(&dest[j],current_alpha,&src[j]);
   }

   free(al_active_ind);
   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorInnerProd        v=x'*y  using indices        complex
 *
 * call seq < MultiInnerProd < i->MultiInnerProd < mv_MultiVectorByMultiVector
 *--------------------------------------------------------------------------*/
int serial_Multi_VectorInnerProd( serial_Multi_Vector *x,
                                   serial_Multi_Vector *y,
                                   int gh, int h, int w, komplex* v)
{
   /* to be reworked! */
   komplex *x_data;
   komplex *y_data;
   int      size;
   int      x_num_active_vectors;
   int      y_num_active_vectors;
   int      i,j,k;
   komplex *y_ptr;
   komplex *x_ptr;
   int * x_active_ind;
   int * y_active_ind;
   komplex current_product;
   komplex temp;
   komplex conj;
   int gap;

   assert (x->size==y->size);

   x_data = (komplex *) x->data;
   y_data = (komplex *) y->data;
   size = x->size;
   x_num_active_vectors = x->num_active_vectors;
   y_num_active_vectors = y->num_active_vectors;

   assert (x_num_active_vectors==h && y_num_active_vectors==w);

   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;

   gap = gh-h;

   for(j=0; j<y_num_active_vectors; j++)
   {
      y_ptr = y_data + y_active_ind[j]*size;

      for (i=0; i<x_num_active_vectors; i++) {

         x_ptr = x_data + x_active_ind[i]*size;
         current_product.real = 0.0;
         current_product.imag = 0.0;

         for(k=0; k<size; k++) {
            conj.real=x_ptr[k].real;
            conj.imag=-x_ptr[k].imag;
            complex_multiply(&temp,&conj,&y_ptr[k]);
            complex_add(&current_product,&current_product,&temp);
         }

         /* fortran column-wise storage for results */
         *v++ = current_product;
      }
      v+=gap;
   }

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorInnerProdDiag                                 complex
 *
 * diag=diagonal(x'*y) using indices
 * mask is index of diag
 *
 * call seq < MultiInnerProdDiag < i->MultiInnerProdDiag < mv_MultiVectorByMultiVectorDiag
 *--------------------------------------------------------------------------*/
int serial_Multi_VectorInnerProdDiag( serial_Multi_Vector *x,
                                       serial_Multi_Vector *y,
                                       int* mask, int n, komplex* diag)
{
/* to be reworked! */
   komplex  *x_data;
   komplex  *y_data;
   int      size;
   int      num_active_vectors;
   int      * x_active_ind;
   int      * y_active_ind;
   komplex  *y_ptr;
   komplex  *x_ptr;
   komplex  current_product;
   komplex  temp;
   komplex  conj;
   int      i, k;
   int      * al_active_ind;
   int      num_active_als;

   assert(x->size==y->size && x->num_active_vectors == y->num_active_vectors);

      /* build list of active indices in alpha */

   al_active_ind = (int *) malloc(sizeof(int)*n);
   num_active_als = 0;

   if (mask!=NULL)
      for (i=0; i<n; i++)
      {
         if (mask[i])
            al_active_ind[num_active_als++]=i;
      }
   else
      for (i=0; i<n; i++)
         al_active_ind[num_active_als++]=i;

   assert (num_active_als==x->num_active_vectors);

   x_data = (komplex *) x->data;
   y_data = (komplex *) y->data;
   size = x->size;
   num_active_vectors = x->num_active_vectors;
   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;

   for (i=0; i<num_active_vectors; i++)
   {
      x_ptr = x_data + x_active_ind[i]*size;
      y_ptr = y_data + y_active_ind[i]*size;
      current_product.real = 0.0;
      current_product.imag = 0.0;

      for(k=0; k<size; k++) {
         conj.real = x_ptr[k].real;
         conj.imag = -x_ptr[k].imag;
         complex_multiply(&temp,&conj,&y_ptr[k]);
         complex_add(&current_product,&current_product,&temp);
      }
      diag[al_active_ind[i]] = current_product;
   }

   free(al_active_ind);
   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorByMatrix      y=x*rVal using indices          complex
 *
 * call seq < MultiVectorByMatrix < i->MultiVecMat < mv_MultiVectorByMatrix
 *--------------------------------------------------------------------------*/
int
 serial_Multi_VectorByMatrix(serial_Multi_Vector *x, int rGHeight, int rHeight,
                              int rWidth, komplex* rVal, serial_Multi_Vector *y)
{
   komplex *x_data;
   komplex *y_data;
   int      size;
   int     * x_active_ind;
   int     * y_active_ind;
   komplex *y_ptr;
   komplex *x_ptr;
   komplex  current_coef;
   komplex  temp;
   int      i,j,k;
   int      gap;

   assert(rHeight>0);
   assert (rHeight==x->num_active_vectors && rWidth==y->num_active_vectors);

   x_data = (komplex *) x->data;
   y_data = (komplex *) y->data;
   size = x->size;
   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;
   gap = rGHeight - rHeight;

   for (j=0; j<rWidth; j++)
   {
      y_ptr = y_data + y_active_ind[j]*size;

      /* ------ set current "y" to first member in a sum ------ */
      x_ptr = x_data + x_active_ind[0]*size;
      current_coef = *rVal++;

      for (k=0; k<size; k++)
         complex_multiply(&y_ptr[k],&current_coef,&x_ptr[k]);

      /* ------ now add all other members of a sum to "y" ----- */
      for (i=1; i<rHeight; i++)
      {
         x_ptr = x_data + x_active_ind[i]*size;
         current_coef = *rVal++;

         for (k=0; k<size; k++) {
            complex_multiply(&temp,&current_coef,&x_ptr[k]);
            complex_add(&y_ptr[k],&y_ptr[k],&temp);
         }
      }

      rVal += gap;
   }

   return 0;
}
/*--------------------------------------------------------------------------
 * serial_Multi_VectorByMulti_Vector     z=x*y   with indices        complex
 *--------------------------------------------------------------------------*/
int
 serial_Multi_VectorByMulti_Vector(serial_Multi_Vector *x,
                                   serial_Multi_Vector *y,
                                   serial_Multi_Vector *z)
{
   komplex *x_data;
   komplex *y_data;
   komplex *z_data;

   int     * x_index;
   int     * y_index;
   int     * z_index;
   komplex * pzc;
   komplex * pyc;
   komplex * pxr;
   komplex * py;

   komplex  temp;
   int      i,j,k;

   assert (x->num_active_vectors == y->size);
   assert (z->size == x->size);
   assert (z->num_active_vectors == y->num_active_vectors);

   x_data = (komplex *) x->data;
   y_data = (komplex *) y->data;
   z_data = (komplex *) z->data;

   x_index = x->active_indices;
   y_index = y->active_indices;
   z_index = z->active_indices;

   for (j=0; j<y->num_active_vectors; j++) {
      pzc = z_data + z->size*z_index[j];
      pyc = y_data + y->size*y_index[j];
      for (i=0;i<x->size;i++) {
         pxr = x_data + i + x->size*x_index[0];
         pzc->real = 0.0;
         pzc->imag = 0.0;
         for (k=0, py=pyc; k<y->size; k++) {
            complex_multiply(&temp, pxr, py);
            complex_add(pzc, pzc, &temp);
            pxr += x->size;
            py++;
         }
         pzc++;
      }
   }

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorPrint                                          complex
 *--------------------------------------------------------------------------*/
int serial_Multi_VectorPrint(serial_Multi_Vector * x,char * tag, int limit)
{
   komplex * p;
   int     * pact;
   int       i, j;
   int     rows,cols;
   printf("======= %s =========\n",tag);
   printf("size %d\n", x->size);
   printf("owns data %d\n",x->owns_data);
   printf("num vectors %d\n",x->num_vectors);
   printf("num active vectors  %d\n",x->num_active_vectors);

   rows = x->size;
   cols = x->num_active_vectors;
   if (limit != 0) {
     if (rows > limit) rows = limit;
     if (cols > limit) cols = limit;
   }

   pact=x->active_indices;
   for (i=0; i<cols; i++, pact++)
      printf("index %d active %d\n", i, *pact);

   for (i=0; i<cols; i++)
   {  p=(komplex *)x->data;
      p = &p[x->active_indices[i]*x->size];
      for (j = 0; j < rows; j++,p++)
         printf("%d %d  %22.16e  %22.16e \n",j,i,p->real,p->imag);
   }

   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorPrintShort                                     complex
 *--------------------------------------------------------------------------*/
int serial_Multi_VectorPrintShort(serial_Multi_Vector * x)
{
   printf("size %d\n", x->size);
   printf("owns data %d\n",x->owns_data);
   printf("num vectors %d\n",x->num_vectors);
   printf("num active vectors  %d\n",x->num_active_vectors);
   return(0);
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorLoad                                          complex
 *--------------------------------------------------------------------------*/
serial_Multi_Vector *
serial_Multi_VectorLoad( char fileName[] ) {

int size, num_vectors, num_elements, j;
serial_Multi_Vector *mvector;

FILE* fp;
komplex* p;
float freal, fimag;

if ( (fp=fopen(fileName,"r") )==NULL) {
  printf("file open failed\n");
  return NULL;
}
fscanf(fp,"%d %d",&size, &num_vectors);
num_elements = size*num_vectors;
mvector = serial_Multi_VectorCreate(size, num_vectors);
serial_Multi_VectorInitialize( mvector);

p = (komplex *)mvector->data;
j = 0;
while ( fscanf(fp,"%e %e",&freal,&fimag) != EOF && j<num_elements ) {
    p->real = freal;
    p->imag = fimag;
    j++;
    p++;
}
fclose(fp);
return mvector;
}
/* ------------------------------------------------------------
   serial_Multi_VectorSymmetrize                        complex
                                     mtx=(mtx+transpose(mtx))/2
   ------------------------------------------------------------ */
void
serial_Multi_VectorSymmetrize( serial_Multi_Vector * mtx ) {

  long i, j, g, h, w, jump;
  komplex* p;
  komplex* q;
  komplex conj_q;

  assert( mtx != NULL );

  g = mtx->size;
  h = mtx->size;
  w = mtx->num_vectors;

  assert( h == w );

  jump = 0;

  for ( j = 0, p = (komplex *) mtx->data; j < w; j++ ) {
    q = p;
    for ( i = j; i < h; i++, p++, q += g ) {
      conj_q.real = q->real;
      conj_q.imag = - (q->imag);
      complex_add(p,p,&conj_q);
      p->real = p->real/2;
      p->imag = p->imag/2;
      q->real = p->real;
      q->imag = -(p->imag);
    }
    p += ++jump;
  }
}
