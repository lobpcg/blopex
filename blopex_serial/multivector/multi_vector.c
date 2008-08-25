/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#include "multi_vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*--------------------------------------------------------------------------
 * serial_Multi_VectorCreate
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
 * serial_Multi_VectorInitialize
 *--------------------------------------------------------------------------*/

int 
serial_Multi_VectorInitialize( serial_Multi_Vector *mvector )
{
   int    ierr = 0, i, size, num_vectors;
   
   size        = serial_Multi_VectorSize(mvector);
   num_vectors = serial_Multi_VectorNumVectors(mvector);
   
   if (NULL==serial_Multi_VectorData(mvector))
      serial_Multi_VectorData(mvector) = 
            (double *) malloc (sizeof(double)*size*num_vectors);
	       
   /* now we create a "mask" of "active" vectors; initially all vectors are active */
   if (NULL==mvector->active_indices)
    {
         mvector->active_indices=malloc(sizeof(int)*num_vectors);

         for (i=0; i<num_vectors; i++)
            mvector->active_indices[i]=i;
            
         mvector->num_active_vectors=num_vectors;
    }

   return ierr;
}


/*--------------------------------------------------------------------------
 * serial_Multi_VectorSetDataOwner
 *--------------------------------------------------------------------------*/

int 
serial_Multi_VectorSetDataOwner( serial_Multi_Vector *mvector, int owns_data )
{
   int    ierr=0;

   serial_Multi_VectorOwnsData(mvector) = owns_data;

   return ierr;
}


/*--------------------------------------------------------------------------
 * serial_Multi_VectorDestroy
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
 * serial_Multi_VectorSetMask
 *-------------------------------------------------------------------------*/
 
 int
 serial_Multi_VectorSetMask(serial_Multi_Vector *mvector, int * mask)
 {
   /* this routine accepts mask in "zeros and ones format, and converts it to the one used in
   the structure "serial_Multi_Vector" */
   int  num_vectors = mvector->num_vectors;
   int i;
   
   
   /* may be it's better to just check if it is not null, and throw an error, if it is? */
   if (mvector->active_indices==NULL)
      mvector->active_indices=malloc(sizeof(int)*num_vectors);
      
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
 * serial_Multi_VectorSetConstantValues
 *--------------------------------------------------------------------------*/

int
serial_Multi_VectorSetConstantValues( serial_Multi_Vector *v,
                                        double             value)
{
   double  *vector_data = serial_Multi_VectorData(v);
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
 * serial_Multi_VectorSetRandomValues
 *
 *     returns vector of values randomly distributed between -1.0 and +1.0
 *--------------------------------------------------------------------------*/

int
serial_Multi_VectorSetRandomValues( serial_Multi_Vector *v, int seed)
{
   double  *vector_data = serial_Multi_VectorData(v);
   int      size        = serial_Multi_VectorSize(v);
   int      i, j, start_offset, end_offset;
           
   srand48(seed);

   for (i = 0; i < v->num_active_vectors; i++)
   {   
      start_offset = v->active_indices[i]*size;
      end_offset = start_offset+size;      

      for (j=start_offset; j < end_offset; j++)
         vector_data[j]= 2.0 * drand48() - 1.0;
   }
    
   return 0;
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorCopy
 * copies data from x to y
 * y should have already been initialized at the same size as x
 *--------------------------------------------------------------------------*/

int
serial_Multi_VectorCopy( serial_Multi_Vector *x, serial_Multi_Vector *y)
{
   double  *x_data; 
   double  *y_data;
   int i;
   int size;
   int num_bytes;
   int num_active_vectors;
   double * dest;
   double * src;
   int * x_active_ind;
   int * y_active_ind;
   
   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);
   
   num_active_vectors = x->num_active_vectors;
   size = x->size;
   num_bytes = size*sizeof(double);
   x_data = x->data;
   y_data = y->data;
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
   

int
serial_Multi_VectorCopyWithoutMask(serial_Multi_Vector *x , serial_Multi_Vector *y)
{
   int byte_count;
   
   assert (x->size == y->size && x->num_vectors == y->num_vectors);
   
   byte_count = sizeof(double) * x->size * x->num_vectors;
   
/* threading not done here since it's not done (reason?) in serial_VectorCopy
      from vector.c */

   memcpy(y->data,x->data,byte_count);
   
   return 0;
}
  

/*--------------------------------------------------------------------------
 * serial_Multi_VectorAxpy
 *--------------------------------------------------------------------------*/

int
serial_Multi_VectorAxpy( double            alpha,
                          serial_Multi_Vector *x,
                          serial_Multi_Vector *y)
{
   double  *x_data; 
   double  *y_data;
   double * src;
   double * dest;
   int * x_active_ind;
   int * y_active_ind;
   int i, j;
   int size;
   int num_active_vectors;

   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);
   
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

      for (j=0; j<size; j++)
         dest[j] += alpha*src[j];
   }
   
   return 0;
}
   

/*--------------------------------------------------------------------------
 * serial_Multi_VectorByDiag: " y(<y_mask>) = alpha(<mask>) .* x(<x_mask>) "
 *--------------------------------------------------------------------------*/

int
serial_Multi_VectorByDiag( serial_Multi_Vector *x,
                             int                *mask, 
                             int                n,
                             double             *alpha,
                             serial_Multi_Vector *y)
{
   double  *x_data;
   double  *y_data;
   int      size;
   int      num_active_vectors;
   int      i,j;
   double  *dest;
   double  *src;
   int * x_active_ind;
   int * y_active_ind;
   int * al_active_ind;
   int num_active_als;
   double current_alpha;

   assert (x->size == y->size && x->num_active_vectors == y->num_active_vectors);
      
   /* build list of active indices in alpha */
   
   al_active_ind = malloc(sizeof(int)*n);
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
      current_alpha=alpha[ al_active_ind[i] ];

      for (j=0; j<size; j++)
         dest[j] = current_alpha*src[j];
   }

   free(al_active_ind);
   return 0; 
}

/*--------------------------------------------------------------------------
 * serial_Multi_VectorInnerProd
 *--------------------------------------------------------------------------*/

int serial_Multi_VectorInnerProd( serial_Multi_Vector *x,
                                  serial_Multi_Vector *y, 
                                  int gh, int h, int w, double* v)
{
   /* to be reworked! */
   double  *x_data;
   double  *y_data; 
   int      size;
   int      x_num_active_vectors;
   int      y_num_active_vectors;
   int      i,j,k;
   double  *y_ptr;
   double  *x_ptr;
   int * x_active_ind;
   int * y_active_ind;
   double current_product;
   int gap;
   
   assert (x->size==y->size);

   x_data = x->data;
   y_data = y->data;
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

      for (i=0; i<x_num_active_vectors; i++)
      {
         x_ptr = x_data + x_active_ind[i]*size;
         current_product = 0.0;

         for(k=0; k<size; k++)
            current_product += x_ptr[k]*y_ptr[k];
         
         /* fortran column-wise storage for results */
            *v++ = current_product;
      }
      v+=gap; 
   }

   return 0;	    
}
   
   
/*--------------------------------------------------------------------------
 * serial_Multi_VectorInnerProdDiag  
 *--------------------------------------------------------------------------*/

int serial_Multi_VectorInnerProdDiag( serial_Multi_Vector *x,
                                      serial_Multi_Vector *y, 
				      int* mask, int n, double* diag)
{
/* to be reworked! */
   double   *x_data;
   double   *y_data;
   int      size;
   int      num_active_vectors;
   int      * x_active_ind;
   int      * y_active_ind;  
   double   *y_ptr;
   double   *x_ptr;
   double   current_product;
   int      i, k;
   int      * al_active_ind;
   int      num_active_als;
        
   assert(x->size==y->size && x->num_active_vectors == y->num_active_vectors);
   
      /* build list of active indices in alpha */
   
   al_active_ind = malloc(sizeof(int)*n);
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
   
   x_data = x->data;
   y_data = y->data;
   size = x->size;
   num_active_vectors = x->num_active_vectors;
   x_active_ind = x->active_indices;
   y_active_ind = y->active_indices;   
   
   for (i=0; i<num_active_vectors; i++)
   {
      x_ptr = x_data + x_active_ind[i]*size;
      y_ptr = y_data + y_active_ind[i]*size;
      current_product = 0.0;
      
      for(k=0; k<size; k++)
            current_product += x_ptr[k]*y_ptr[k];
      
      diag[al_active_ind[i]] = current_product;      
   }
    
   free(al_active_ind);
   return 0;	    
}
   

int
 serial_Multi_VectorByMatrix(serial_Multi_Vector *x, int rGHeight, int rHeight, 
                                 int rWidth, double* rVal, serial_Multi_Vector *y)
{
   double   *x_data;
   double   *y_data;
   int      size;
   int      * x_active_ind;
   int      * y_active_ind;  
   double   *y_ptr;
   double   *x_ptr;
   double   current_coef;
   int      i,j,k;
   int      gap;
      
   assert(rHeight>0);
   assert (rHeight==x->num_active_vectors && rWidth==y->num_active_vectors);
   
   x_data = x->data;
   y_data = y->data;
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
         y_ptr[k] = current_coef * x_ptr[k];
         
      /* ------ now add all other members of a sum to "y" ----- */
      for (i=1; i<rHeight; i++)
      {
         x_ptr = x_data + x_active_ind[i]*size;
         current_coef = *rVal++;

         for (k=0; k<size; k++)
            y_ptr[k] += current_coef * x_ptr[k];
      }

      rVal += gap;
   }

   return 0;
}

int serial_Multi_VectorPrint(serial_Multi_Vector * x, const char *base_name)
{
   FILE     *fp;
   int      i, j;
   char     fullName[128];
   
   for (i=0; i<x->num_vectors; i++)
   {
      sprintf( fullName, "%s.%d", base_name, i );
      fp = fopen(fullName, "w");
      assert (fp!=NULL);
      
      fprintf(fp, "%d\n", x->size);
      for (j = 0; j < x->size; j++)
         fprintf(fp, "%e\n", x->data[j + i*x->size] );
      
      fclose(fp);
   }

   return 0;
}
