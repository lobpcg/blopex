/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#ifndef serial_MULTI_VECTOR_HEADER
#define serial_MULTI_VECTOR_HEADER

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * serial_Multi_Vector
 *--------------------------------------------------------------------------*/

typedef struct
{
   double  *data;
   int      size;
   int      owns_data;
   int      num_vectors;  /* the above "size" is size of one vector */

   int      num_active_vectors;
   int     *active_indices;   /* indices of active vectors; 0-based notation */

} serial_Multi_Vector ;

/*--------------------------------------------------------------------------
 * Accessor functions for the Multi_Vector structure
 *--------------------------------------------------------------------------*/

#define serial_Multi_VectorData(vector)      ((vector) -> data)
#define serial_Multi_VectorSize(vector)      ((vector) -> size)
#define serial_Multi_VectorOwnsData(vector)  ((vector) -> owns_data)
#define serial_Multi_VectorNumVectors(vector) ((vector) -> num_vectors)

serial_Multi_Vector * serial_Multi_VectorCreate( int size, int num_vectors  );
serial_Multi_Vector *serial_Multi_VectorRead( char *file_name );

int serial_Multi_VectorDestroy( serial_Multi_Vector *vector );
int serial_Multi_VectorInitialize( serial_Multi_Vector *vector );
int serial_Multi_VectorSetDataOwner(serial_Multi_Vector *vector , int owns_data);
/*
int serial_Multi_VectorPrint( serial_Multi_Vector *vector, const char *file_name );
*/
int serial_Multi_VectorPrint( serial_Multi_Vector *vector, char * tag, int limit);

int serial_Multi_VectorSetConstantValues(serial_Multi_Vector *v,double value);
int serial_Multi_VectorSetRandomValues(serial_Multi_Vector *v , int seed);
int serial_Multi_VectorCopy( serial_Multi_Vector *x , serial_Multi_Vector *y);
int serial_Multi_VectorScale( double alpha , serial_Multi_Vector *y, int *mask  );
int serial_Multi_VectorAxpy( double alpha , serial_Multi_Vector *x , serial_Multi_Vector *y);
int serial_Multi_VectorInnerProd( serial_Multi_Vector *x,
                                  serial_Multi_Vector *y,
                                  int gh, int h, int w, double* v);
int serial_Multi_VectorMultiScale( double *alpha, serial_Multi_Vector *v, int *mask );

int serial_Multi_VectorByDiag( serial_Multi_Vector *x,
                                 int                *mask,
                                 int                n,
                                 double             *alpha,
                                 serial_Multi_Vector *y);

int serial_Multi_VectorInnerProdDiag( serial_Multi_Vector *x,
                                      serial_Multi_Vector *y,
                      int* mask, int n, double* diag);

int
 serial_Multi_VectorSetMask(serial_Multi_Vector *mvector, int * mask);
int
 serial_Multi_VectorCopyWithoutMask(serial_Multi_Vector *x , serial_Multi_Vector *y);
int
 serial_Multi_VectorByMatrix(serial_Multi_Vector *x, int rGHeight, int rHeight,
                                 int rWidth, double* rVal, serial_Multi_Vector *y);
int
 serial_Multi_VectorByMulti_Vector(serial_Multi_Vector *x,
                                   serial_Multi_Vector *y,
                                   serial_Multi_Vector *z);
#ifdef __cplusplus
}
#endif

#endif /* serial_MULTI_VECTOR_HEADER */
