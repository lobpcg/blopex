/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* @@@ BLOPEX (version 1.1) LGPL Version 2.1 or above.See www.gnu.org. */
/* @@@ Copyright 2010 BLOPEX team http://code.google.com/p/blopex/     */
/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#include "matmultivec.h"
#include "multi_vector.h"
#include <assert.h>

BlopexInt serial_Lapl1DMatMultiVec( serial_Multi_Vector * x,
                            serial_Multi_Vector * y  )
{
   /* this function implements 1D Laplacian operator, Dirichlet BC */

   double  *x_data; 
   double  *y_data;
   double * src;
   double * dest;
   BlopexInt * x_active_ind;
   BlopexInt * y_active_ind;
   BlopexInt i, j;
   BlopexInt size;
   BlopexInt num_active_vectors;

   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);
   assert (x->size>1);
   
   x_data = x->data;
   y_data = y->data;
   size = x->size;
   num_active_vectors = x->num_active_vectors;
   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;
   
   for(i=0; i<num_active_vectors; i++)
   {
      src = x_data + x_active_ind[i]*size;
      dest = y_data + y_active_ind[i]*size;

      dest[0] = 2*src[0] - src[1];
      for (j=1; j<size-1; j++)
         dest[j] = -src[j-1] + 2*src[j] - src[j+1];
      dest[size-1] = -src[size-2] + 2*src[size-1];
   }
   
   return 0;
}                            
